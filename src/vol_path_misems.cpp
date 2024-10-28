#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

#define MINIMUM_BOUNCES 3
#define MIN(a, b) (((a)>(b))?(b):(a))

class VolPathMisEms : public Integrator {
public:
    VolPathMisEms(const PropertyList &props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Iterative algorithm */

        int bounces = 0;
        float etaProd = 1.0f;
        float continueProbability = 1.0f;
        bool continueBouncing = true;
        bool isSpecularBounce = true;
        Color3f emittedAndReflected = Color3f(0.0f, 0.0f, 0.0f);
        Color3f throughput = Color3f(1.0f, 1.0f, 1.0f);
        Ray3f* _ray = new Ray3f(ray);
        Medium* nextMedium = scene->getMedium();

        Intersection nextIts;

        /* Find the surface that is visible in the requested direction */
        bool nextHitSurface = scene->rayIntersect(*_ray, nextIts);

        if (!nextHitSurface) {
            // If no medium, usual termination
            if (nextMedium == nullptr) {
                // May have hit an environment map!
                if (isSpecularBounce && scene->getEnvironmentMap()) {
                    Color3f emitted = scene->getEnvironmentMap()->eval(_ray->d);
                    emittedAndReflected += (throughput * emitted / continueProbability);
                }
                    
                delete _ray;
                return emittedAndReflected;
            }
        } else {
            _ray->maxt = nextIts.t;
        }


        // If there is a medium, let the medium scatter
        Point3f nextInteractionPoint;
        bool nextSampledMedium = false;
        Color3f nextFr = Color3f(1.0f);
        Point3f nextNewPoint;
        if (nextMedium != nullptr) {
            nextFr = nextMedium->sample(*_ray, *sampler, nextInteractionPoint, nextSampledMedium);
            if (nextFr.cwiseEqual(Color3f(0.0f)).all()) {
                delete _ray;
                return emittedAndReflected;
            }
        }
        if (nextSampledMedium) {
            nextNewPoint = nextInteractionPoint;
        } else {
            nextNewPoint = nextIts.p;
        }

        while (continueBouncing) {

            /* Find the surface that is visible in the requested direction */
            Intersection its(nextIts);
            bool hitSurface = nextHitSurface;
            Medium *medium = nextMedium;

            // If there is a medium, let the medium scatter
            bool sampledMedium = nextSampledMedium;
            Color3f fr = nextFr;
            Point3f newPoint = nextNewPoint;

            /* Sample direct illumination (only if in a medium or a non-dirac surface)*/
            if (!scene->getEmitters().empty() && (sampledMedium || (hitSurface && its.mesh->getBSDF()->isDiffuse()))) {
                Point2f sample = sampler->next2D();
                uint32_t emitterIndex = floor(sampler->next1D() * scene->getEmitters().size());

                Point3f samplePosition;
                Normal3f sampleNormal;

                float epsilon;
                Color3f le = scene->getEmitters()[emitterIndex]->sample(newPoint, sample, samplePosition, sampleNormal, epsilon);
                if (!(sampleNormal.dot(newPoint - samplePosition) < 0) && (sampledMedium || !(Frame::cosTheta(its.toLocal(samplePosition - newPoint)) <= 0))) {
                    //Instead of checking just one shadow ray, I cycle through ray intersections, because I don't care about null surfaces that only change the medium
                    // While doing so, I accumulate the transmittance of the media

                    Color3f transmittance = Color3f(1.0f);
                    Point3f currentShadowPoint = newPoint;
                    Medium *shadowMedium = medium;
                    // Start from the current medium, but if we hit a surface, start from the medium of the surface instead
                    if (!sampledMedium) {
                        if (its.mesh->hasMedia()) {
                            if (its.shFrame.cosTheta((samplePosition - newPoint).normalized()) > 0.0f) { //Outgoing direction agrees with the surface normal, get the external medium
                                shadowMedium = its.mesh->getMediumExternal();
                            }   
                            else { //Get the internal medium
                                shadowMedium = its.mesh->getMediumInternal();
                            }
                        }
                    }

                    while (true) {
                        Ray3f shadowRay = Ray3f(currentShadowPoint, (samplePosition - currentShadowPoint));
                        shadowRay.mint = epsilon;
                        shadowRay.maxt = 1 - epsilon;
                        Intersection shadowIts;
                        bool shadowHitSurface = scene->rayIntersect(shadowRay, shadowIts);
                        
                        if (shadowHitSurface) {
                            if (shadowIts.mesh->getBSDF()->isNull()) {
                                // We've encountered a null surface: update the point, medium, transmittance, and keep going
                                shadowRay.maxt = shadowIts.t;
                                if (shadowMedium != nullptr)
                                    transmittance *= shadowMedium->transmittance(shadowRay, *sampler);
                                if (shadowIts.mesh->hasMedia()) {
                                    if (shadowIts.shFrame.cosTheta(shadowIts.toLocal((samplePosition - newPoint).normalized())) > 0.0f) { 
                                        //Outgoing direction agrees with the surface normal, get the external medium
                                        shadowMedium = shadowIts.mesh->getMediumExternal();
                                    }   
                                    else { //Get the internal medium
                                        shadowMedium = shadowIts.mesh->getMediumInternal();
                                    }
                                }
                                currentShadowPoint = shadowIts.p;
                            } else 
                                break; // We've encountered a solid object, no direct light from this light sample
                        } else {
                            // We didn't encounter an object, but we should compute transmittance anyway
                            if (shadowMedium != nullptr)
                                transmittance *= shadowMedium->transmittance(shadowRay, *sampler);

                            // If at some point there is no intersection, compute as usual and account for the media transmittance in shadowFr
                            if (sampledMedium) {
                                Color3f shadowFr = medium->getPhaseFunction()->pdf(-_ray->d.normalized(), (samplePosition - newPoint).normalized());

                                Vector3f distance = samplePosition - newPoint;

                                float p_light = scene->getEmitters()[emitterIndex]->pdf(newPoint, samplePosition) / (sampleNormal.normalized().dot((newPoint - samplePosition).normalized())) / scene->getEmitters().size();
                                float p_phase = 0.0f;
                                if (scene->getEmitters()[emitterIndex]->canBeHit())
                                    p_phase = medium->getPhaseFunction()->pdf(-_ray->d.normalized(), (samplePosition - newPoint).normalized());
                                if (scene->getEmitters()[emitterIndex]->dependsOnDistance())
                                    p_light *= distance.squaredNorm();
                                // float mis_weight = p_light / ((p_light + p_phase));
                                float mis_weight = 1.0f;

                                float g = abs(sampleNormal.dot((samplePosition - newPoint).normalized())) / (newPoint - samplePosition).squaredNorm();

                                // if (p_light >= 0 && p_phase >= 0) {
                                    emittedAndReflected += shadowFr.cwiseProduct(le).cwiseProduct(throughput) * g * fr * transmittance * scene->getEmitters().size() / continueProbability * mis_weight;
                                // }
                            } else {
                                BSDFQueryRecord query = BSDFQueryRecord(its.toLocal(-_ray->d.normalized()), its.toLocal((samplePosition - newPoint).normalized()), ESolidAngle, &its);

                                Color3f shadowFr = its.mesh->getBSDF()->eval(query);

                                Vector3f distance = samplePosition - newPoint;

                                float p_light = scene->getEmitters()[emitterIndex]->pdf(newPoint, samplePosition) / (sampleNormal.normalized().dot((newPoint - samplePosition).normalized())) / scene->getEmitters().size();
                                float p_brdf = 0.0f;
                                if (scene->getEmitters()[emitterIndex]->canBeHit())
                                    p_brdf = its.mesh->getBSDF()->pdf(query);
                                if (scene->getEmitters()[emitterIndex]->dependsOnDistance())
                                    p_light *= distance.squaredNorm();
                                // float mis_weight = p_light / ((p_light + p_brdf));
                                float mis_weight = 1.0f;


                                float g = abs(its.shFrame.cosTheta(its.toLocal((samplePosition - newPoint).normalized())) * sampleNormal.dot((samplePosition - newPoint).normalized())) / (newPoint - samplePosition).squaredNorm();
                                
                                // if (p_light >= 0 && p_brdf >= 0) {
                                    emittedAndReflected += shadowFr.cwiseProduct(le).cwiseProduct(throughput) * g * fr * transmittance * scene->getEmitters().size() / continueProbability * mis_weight;
                                // }
                            }
                            break;
                        }
                    }
                }
            }
            


            Vector3f newDirection;
            float eta = 1.0f;
            //If the medium scattered the light
            if (sampledMedium) {
                // float phaseFr = medium->getPhaseFunction()->sample(-_ray->d.normalized(), &newDirection, sampler->next2D());
                medium->getPhaseFunction()->sample(-_ray->d.normalized(), &newDirection, sampler->next2D());
                //No modification to the fr value as phaseFr/phasePDF = 1
            }
            //Otherwise proceed like in the usual path_mis, but considering the transmittance of the medium as well
            else {

                if (!hitSurface) {
                    // May have hit an environment map!
                    if (isSpecularBounce && scene->getEnvironmentMap()) {
                        Color3f emitted = scene->getEnvironmentMap()->eval(_ray->d);
                        emittedAndReflected += (throughput * emitted * fr / continueProbability);
                    }

                    delete _ray;
                    return emittedAndReflected;
                }


                if (isSpecularBounce && its.mesh->isEmitter()) {
                    Color3f emitted = its.mesh->getEmitter()->eval(its.toLocal(-_ray->d).normalized(), &its);
                    emittedAndReflected += (throughput * emitted * fr / continueProbability);
                }


                /* Sampling a new direction */
                BSDFQueryRecord query = BSDFQueryRecord(its.toLocal(-_ray->d).normalized(), sampler, &its);
                Color3f bsdfFr = its.mesh->getBSDF()->sample(query, sampler->next2D());
                
                if (bsdfFr.cwiseEqual(Color3f(0.0f)).all()) {
                    delete _ray;
                    return emittedAndReflected;
                }

                // We've hit an object, update the medium if the object has any media
                if (its.mesh->hasMedia()) {
                    if (its.shFrame.cosTheta(query.wo) > 0.0f) { //Outgoing direction agrees with the surface normal, get the external medium
                        nextMedium = its.mesh->getMediumExternal();
                    }   
                    else { //Get the internal medium
                        nextMedium = its.mesh->getMediumInternal();
                    }
                }

                fr *= bsdfFr;
                newDirection = its.toWorld(query.wo);
                eta = query.eta;

            }


            /* Checking if the new direction has a light source, and storing the intersection anyway for the next iteration */
            Ray3f nextRay = Ray3f(newPoint, newDirection);
            nextRay.mint = Epsilon;
            nextHitSurface = scene->rayIntersect(nextRay, nextIts);
            if (nextHitSurface)
                nextRay.maxt = nextIts.t;

            nextSampledMedium = false;
            nextFr = Color3f(1.0f);
            if (nextMedium != nullptr) {
                nextFr = nextMedium->sample(nextRay, *sampler, nextInteractionPoint, nextSampledMedium);
                if (nextFr.cwiseEqual(Color3f(0.0f)).all()) {
                    delete _ray;
                    return emittedAndReflected;
                }
            }
            if (nextSampledMedium) {
                nextNewPoint = nextInteractionPoint;
            } else {
                nextNewPoint = nextIts.p;
            }

            bool nextIsSpecularBounce;
            if (sampledMedium)
                nextIsSpecularBounce = false; //Media, for now, cannot be Dirac's deltas.
            else
                if (its.mesh->getBSDF()->isNull()) {
                    nextIsSpecularBounce = isSpecularBounce;
                } else
                    nextIsSpecularBounce = !its.mesh->getBSDF()->isDiffuse();

            if (nextHitSurface) {
                if (!nextSampledMedium) {
                    if (nextIts.mesh->isEmitter()) {
                        Color3f emitted = nextIts.mesh->getEmitter()->eval(nextIts.toLocal(-nextRay.d).normalized(), &nextIts);

                        Vector3f distance = nextIts.p - newPoint;

                        float p_light = nextIts.mesh->getEmitter()->pdf(newPoint, nextIts.toLocal(-nextRay.d)) * distance.squaredNorm() / Frame::cosTheta(nextIts.toLocal((its.p - nextIts.p).normalized())) / scene->getEmitters().size();
                        
                        if (sampledMedium) {
                            float p_phase = medium->getPhaseFunction()->pdf(-ray.d.normalized(), newDirection);
                            // float mis_weight = p_phase / (p_light + p_phase);
                            float mis_weight = 0.0f;
                            // if (p_light >= 0 && p_phase > 0)
                                emittedAndReflected += (throughput * emitted * fr * nextFr / continueProbability) * mis_weight;

                        } else {
                            BSDFQueryRecord query = BSDFQueryRecord(its.toLocal(-_ray->d).normalized(), newDirection, ESolidAngle, sampler, &its);
                            float p_brdf = its.mesh->getBSDF()->pdf(query);
                            // float mis_weight = p_brdf / ((p_light + p_brdf));
                            float mis_weight = 0.0f;
                            // if (p_light >= 0 && p_brdf > 0)
                            if (!nextIsSpecularBounce)
                                emittedAndReflected += (throughput * emitted * fr * nextFr / continueProbability) * mis_weight;
                        }
                    
                    }
                }
                
            } else {
                if (!nextSampledMedium) {
                    if (scene->getEnvironmentMap()) {
                        Color3f emitted = scene->getEnvironmentMap()->eval(newDirection);
                        float p_light = scene->getEnvironmentMap()->pdf(newPoint, newPoint + newDirection) / scene->getEmitters().size();

                        if (sampledMedium) {
                            float p_phase = medium->getPhaseFunction()->pdf(-ray.d.normalized(), newDirection);
                            // float mis_weight = p_phase / (p_light + p_phase);
                            float mis_weight = 0.0f;
                            // if (p_light >= 0 && p_phase > 0)
                                emitted *= mis_weight;
                            // else
                                // emitted = Color3f(0.0f);
                            emittedAndReflected += (throughput * emitted * fr * nextFr / continueProbability) * mis_weight;
                        } else {
                            if (its.mesh->getBSDF()->isDiffuse()) {
                                BSDFQueryRecord query = BSDFQueryRecord(its.toLocal(-_ray->d).normalized(), newDirection, ESolidAngle, sampler, &its);
                                float p_brdf = its.mesh->getBSDF()->pdf(query);
                                // float mis_weight = p_brdf / ((p_light + p_brdf));
                                float mis_weight = 0.0f;
                                // if (p_light >= 0 && p_brdf > 0)
                                    emitted *= mis_weight;
                                // else
                                    // emitted = Color3f(0.0f);
                            }
                            emittedAndReflected += (throughput * emitted * fr * nextFr / continueProbability);
                        }
                    }

                    delete _ray;
                    return emittedAndReflected;
                }
                
            }


            delete _ray;
            _ray = new Ray3f(newPoint, newDirection);
            // _ray->mint = Epsilon;

            if (bounces == 0)
                throughput = fr;
            else 
                throughput = throughput * fr / continueProbability;

            etaProd *= eta;
            continueProbability = bounces < MINIMUM_BOUNCES ? 1.0f : MIN(throughput.maxCoeff() * etaProd * etaProd, 0.99f);
            continueBouncing = (bounces < MINIMUM_BOUNCES) || (sampler->next1D() < continueProbability);

            isSpecularBounce = nextIsSpecularBounce;

            bounces += 1;
        }

        delete _ray;
        return emittedAndReflected;
    }

    std::string toString() const {
        return "VolPathMisEms[]";
    }

};

NORI_REGISTER_CLASS(VolPathMisEms, "vol_path_misems");
NORI_NAMESPACE_END