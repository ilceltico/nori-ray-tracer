#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

#define MINIMUM_BOUNCES 3
#define MIN(a, b) (((a)>(b))?(b):(a))

class PathMis : public Integrator {
public:
    PathMis(const PropertyList &props) {
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

        Intersection nextIts;

        if (!scene->rayIntersect(*_ray, nextIts)) {
            // May have hit an environment map!
            if (isSpecularBounce && scene->getEnvironmentMap()) {
                Color3f emitted = scene->getEnvironmentMap()->eval(_ray->d);
                emittedAndReflected += (throughput * emitted / continueProbability);
            }
                
            delete _ray;
            return emittedAndReflected;
        }

        while (continueBouncing) {

            /* Find the surface that is visible in the requested direction */
            Intersection its(nextIts);

            /* Sample direct illumination */
            if (!scene->getEmitters().empty()) {
                Point2f sample = sampler->next2D();
                uint32_t emitterIndex = floor(sampler->next1D() * scene->getEmitters().size());

                Point3f samplePosition;
                Normal3f sampleNormal;

                float epsilon;
                Color3f le = scene->getEmitters()[emitterIndex]->sample(its.p, sample, samplePosition, sampleNormal, epsilon);
                if (! (sampleNormal.dot(its.p - samplePosition) <= 0) && ! (Frame::cosTheta(its.toLocal(samplePosition - its.p)) <= 0) ) {
                    Ray3f shadowRay = Ray3f(its.p, (samplePosition - its.p));
                    shadowRay.mint = epsilon;
                    shadowRay.maxt = 1 - epsilon;
                    if (!scene->rayIntersect(shadowRay)) {

                        BSDFQueryRecord query = BSDFQueryRecord(its.toLocal(-_ray->d).normalized(), its.toLocal((samplePosition - its.p).normalized()), ESolidAngle, sampler, &its);

                        Color3f fr = its.mesh->getBSDF()->eval(query);
                        
                        Vector3f distance = samplePosition - its.p;

                        float p_light = scene->getEmitters()[emitterIndex]->pdf(its.p, samplePosition) / (sampleNormal.normalized().dot((its.p - samplePosition).normalized())) / scene->getEmitters().size();
                        float p_brdf = 0.0f;
                        if (scene->getEmitters()[emitterIndex]->canBeHit())
                            p_brdf = its.mesh->getBSDF()->pdf(query);
                        if (scene->getEmitters()[emitterIndex]->dependsOnDistance())
                            p_light *= distance.squaredNorm();
                        float mis_weight = p_light / ((p_light + p_brdf));
 
                        float g = abs(its.shFrame.cosTheta(its.toLocal((samplePosition - its.p).normalized())) * sampleNormal.dot((its.p - samplePosition).normalized())) / (its.p - samplePosition).squaredNorm();
 
                        if (p_light >= 0 && p_brdf >= 0) {
                            emittedAndReflected += fr.cwiseProduct(le).cwiseProduct(throughput) * g * scene->getEmitters().size() / continueProbability * mis_weight;
                        }
                    }
                }
            }
            

            if (isSpecularBounce && its.mesh->isEmitter()) {
                Color3f emitted = its.mesh->getEmitter()->eval(its.toLocal(-_ray->d).normalized(), &its);
                emittedAndReflected += (throughput * emitted / continueProbability);
            }


            /* Sampling a new direction */
            BSDFQueryRecord query = BSDFQueryRecord(its.toLocal(-_ray->d).normalized(), sampler, &its);
            Color3f fr = its.mesh->getBSDF()->sample(query, sampler->next2D());
            
            if (fr.cwiseEqual(Color3f(0.0f)).all()) {
                delete _ray;
                return emittedAndReflected;
            }

            /* Checking if it's a light source, and storing the intersection anyway for the next iteration */
            if (scene->rayIntersect(Ray3f(its.p, its.toWorld(query.wo)), nextIts)) {
                if (nextIts.mesh->isEmitter()) {
                    Color3f emitted = nextIts.mesh->getEmitter()->eval(nextIts.toLocal(-its.toWorld(query.wo)).normalized(), &nextIts);

                    Vector3f distance = nextIts.p - its.p;

                    float p_light = nextIts.mesh->getEmitter()->pdf(its.p, nextIts.toLocal(-its.toWorld(query.wo))) * distance.squaredNorm() / Frame::cosTheta(nextIts.toLocal((its.p - nextIts.p).normalized())) / scene->getEmitters().size();
                    float p_brdf = its.mesh->getBSDF()->pdf(query);
                    float mis_weight = p_brdf / ((p_light + p_brdf));
                    if (p_light >= 0 && p_brdf > 0)
                        emittedAndReflected += (throughput * emitted * fr / continueProbability) * mis_weight;
                }
            } else {
                if (scene->getEnvironmentMap()) {
                    Color3f emitted = scene->getEnvironmentMap()->eval(its.toWorld(query.wo));
                    if (its.mesh->getBSDF()->isDiffuse()) {
                        float p_light = scene->getEnvironmentMap()->pdf(its.p, its.p + its.toWorld(query.wo)) / scene->getEmitters().size();
                        float p_brdf = its.mesh->getBSDF()->pdf(query);
                        float mis_weight = p_brdf / ((p_light + p_brdf));
                        if (p_light >= 0 && p_brdf > 0)
                            emitted *= mis_weight;
                        else
                            emitted = Color3f(0.0f);
                    }
                    emittedAndReflected += (throughput * emitted * fr / continueProbability);
                }

                delete _ray;
                return emittedAndReflected;
            }


            delete _ray;
            _ray = new Ray3f(its.p, its.toWorld(query.wo));
            _ray->mint = Epsilon;

            if (bounces == 0)
                throughput = fr;
            else 
                throughput = throughput * fr / continueProbability;

            etaProd *= query.eta;
            continueProbability = bounces < MINIMUM_BOUNCES ? 1.0f : MIN(throughput.maxCoeff() * etaProd * etaProd, 0.99f);
            continueBouncing = (bounces < MINIMUM_BOUNCES) || (sampler->next1D() < continueProbability);

            isSpecularBounce = !its.mesh->getBSDF()->isDiffuse();

            bounces += 1;
        }

        delete _ray;
        return emittedAndReflected;
    }

    std::string toString() const {
        return "PathMis[]";
    }

};

NORI_REGISTER_CLASS(PathMis, "path_mis");
NORI_NAMESPACE_END