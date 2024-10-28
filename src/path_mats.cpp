#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

#define MINIMUM_BOUNCES 3
#define MIN(a, b) (((a)>(b))?(b):(a))

class PathMats : public Integrator {
public:
    PathMats(const PropertyList &props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Iterative algorithm */

        int bounces = 0;
        float etaProd = 1.0f;
        float continueProbability = 1.0f;
        bool continueBouncing = true;
        Color3f emittedAndReflected = Color3f(0.0f, 0.0f, 0.0f);
        Color3f throughput = Color3f(1.0f, 1.0f, 1.0f);
        Ray3f* _ray = new Ray3f(ray);
        while (continueBouncing) {

            /* Find the surface that is visible in the requested direction */
            Intersection its;
            if (!scene->rayIntersect(*_ray, its)) {
                // May have hit an environment map!
                if (scene->getEnvironmentMap()) {
                    Color3f emitted = scene->getEnvironmentMap()->eval(_ray->d);
                    emittedAndReflected += (throughput * emitted / continueProbability);
                }

                delete _ray;
                return emittedAndReflected;
            }

            if (its.mesh->isEmitter()) {
                Color3f emitted = its.mesh->getEmitter()->eval(its.toLocal(-_ray->d.normalized()), &its);
                emittedAndReflected += (throughput * emitted / continueProbability);
            }

            BSDFQueryRecord query = BSDFQueryRecord(its.toLocal(-_ray->d.normalized()), sampler, &its);
            Color3f fr = its.mesh->getBSDF()->sample(query, sampler->next2D());
            if (fr.cwiseEqual(Color3f(0.0f)).all()) {
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

            bounces += 1;
        }

        delete _ray;
        return emittedAndReflected;
    }

    std::string toString() const {
        return "PathMats[]";
    }

};

NORI_REGISTER_CLASS(PathMats, "path_mats");
NORI_NAMESPACE_END