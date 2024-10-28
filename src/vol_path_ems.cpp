#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

#define MINIMUM_BOUNCES 3
#define MIN(a, b) (((a)>(b))?(b):(a))

class VolPathEms : public Integrator {
public:
    VolPathEms(const PropertyList &props) {
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
        Medium* medium = scene->getMedium();
        while (continueBouncing) {

            /* Find the next surface (it can also be null, as it can mean a medium change) */
            Intersection its;
            bool hitSurface = scene->rayIntersect(*_ray, its);
            if (!hitSurface) {
                // If no medium, usual termination
                if (medium == nullptr) {
                    // May have hit an environment map!
                    if (isSpecularBounce && scene->getEnvironmentMap()) {
                        Color3f emitted = scene->getEnvironmentMap()->eval(_ray->d);
                        emittedAndReflected += (throughput * emitted / continueProbability);
                    }

                    delete _ray;
                    return emittedAndReflected;
                }
            } else {
                _ray->maxt = its.t;
            }

            // If there is a medium, let the medium scatter
            Point3f interactionPoint;
            bool sampledMedium = false;
            Color3f fr = Color3f(1.0f);
            Point3f newPoint;
            if (medium != nullptr) {
                fr = medium->sample(*_ray, *sampler, interactionPoint, sampledMedium);
                if (fr.cwiseEqual(Color3f(0.0f)).all()) {
                    delete _ray;
                    return emittedAndReflected;
                }
            }
            if (sampledMedium) {
                newPoint = interactionPoint;
            } else {
                newPoint = its.p;
            }


            /* Sample direct illumination (only if in a medium or a non-dirac surface)*/
            if (!scene->getEmitters().empty() && (sampledMedium || (hitSurface && its.mesh->getBSDF()->isDiffuse()))) {
                Point2f sample = sampler->next2D();
                uint32_t emitterIndex = floor(sampler->next1D() * scene->getEmitters().size());

                Point3f samplePosition;
                Normal3f sampleNormal;

                float epsilon;
                Color3f le = scene->getEmitters()[emitterIndex]->sample(newPoint, sample, samplePosition, sampleNormal, epsilon);
                if (! (sampleNormal.dot(newPoint - samplePosition) < 0)) {
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
                        // epsilon = Epsilon;
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
                                float g = abs(sampleNormal.dot((samplePosition - newPoint).normalized())) / (newPoint - samplePosition).squaredNorm();
                                emittedAndReflected += shadowFr.cwiseProduct(le).cwiseProduct(throughput) * g * fr * transmittance * scene->getEmitters().size() / continueProbability;
                            } else {
                                BSDFQueryRecord query = BSDFQueryRecord(its.toLocal(-_ray->d.normalized()), its.toLocal((samplePosition - newPoint).normalized()), ESolidAngle, &its);

                                Color3f shadowFr = its.mesh->getBSDF()->eval(query);
                                float g = abs(its.shFrame.cosTheta(its.toLocal((samplePosition - newPoint).normalized())) * sampleNormal.dot((samplePosition - newPoint).normalized())) / (newPoint - samplePosition).squaredNorm();
                                emittedAndReflected += shadowFr.cwiseProduct(le).cwiseProduct(throughput) * g * fr * transmittance * scene->getEmitters().size() / continueProbability;
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
            //Otherwise proceed like in the usual path_ems, but considering the transmittance of the medium as well
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
                    Color3f emitted = its.mesh->getEmitter()->eval(its.toLocal(-_ray->d.normalized()), &its);
                    emittedAndReflected += (throughput * emitted * fr / continueProbability);
                }

                /* Indirect illumination */

                BSDFQueryRecord query = BSDFQueryRecord(its.toLocal(-_ray->d.normalized()), sampler, &its);
                Color3f bsdfFr = its.mesh->getBSDF()->sample(query, sampler->next2D());
                if (bsdfFr.cwiseEqual(Color3f(0.0f)).all()) {
                    delete _ray;
                    return emittedAndReflected;
                }

                // We've hit an object, update the medium if the object has any media
                if (its.mesh->hasMedia()) {
                    if (its.shFrame.cosTheta(query.wo) > 0.0f) { //Outgoing direction agrees with the surface normal, get the external medium
                        medium = its.mesh->getMediumExternal();
                    }   
                    else { //Get the internal medium
                        medium = its.mesh->getMediumInternal();
                    }
                }

                fr *= bsdfFr;
                newDirection = its.toWorld(query.wo);
                eta = query.eta;
            }
            

            delete _ray;
            _ray = new Ray3f(newPoint, newDirection);
            _ray->mint = Epsilon;

            if (bounces == 0)
                throughput = fr;
            else 
                throughput = throughput * fr / continueProbability;

            etaProd *= eta;
            continueProbability = bounces < MINIMUM_BOUNCES ? 1.0f : MIN(throughput.maxCoeff() * etaProd * etaProd, 0.99f);
            continueBouncing = (bounces < MINIMUM_BOUNCES) || (sampler->next1D() < continueProbability);

            if (sampledMedium)
                isSpecularBounce = false; //Media, for now, cannot be Dirac's deltas.
            else
                if (its.mesh->getBSDF()->isNull()) {
                    //Null BSDF does not modify anything
                } else
                    isSpecularBounce = !its.mesh->getBSDF()->isDiffuse();

            bounces += 1;
        }

        delete _ray;
        return emittedAndReflected;
    }

    std::string toString() const {
        return "VolPathEms[]";
    }

};

NORI_REGISTER_CLASS(VolPathEms, "vol_path_ems");
NORI_NAMESPACE_END