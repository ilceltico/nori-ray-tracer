#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

#define MINIMUM_BOUNCES 3
#define MIN(a, b) (((a)>(b))?(b):(a))

class VolPathMats : public Integrator {
public:
    VolPathMats(const PropertyList &props) {
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
        Medium* medium = scene->getMedium();
        while (continueBouncing) {

            /* Find the next surface (it can also be null, as it can mean a medium change) */
            Intersection its;
            bool hitSurface = scene->rayIntersect(*_ray, its);
            if (!hitSurface) {
                // If no medium, usual termination
                if (medium == nullptr) {
                    // May have hit an environment map!
                    if (scene->getEnvironmentMap()) {
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
            if (medium != nullptr) {
                fr = medium->sample(*_ray, *sampler, interactionPoint, sampledMedium);
                if (fr.cwiseEqual(Color3f(0.0f)).all()) {
                    delete _ray;
                    return emittedAndReflected;
                }
            }
            

            Point3f newPoint;
            Vector3f newDirection;
            float eta = 1.0f;
            //If the medium scattered the light
            if (sampledMedium) {
                newPoint = interactionPoint;
                // float phaseFr = medium->getPhaseFunction()->sample(-_ray->d.normalized(), &newDirection, sampler->next2D());
                medium->getPhaseFunction()->sample(-_ray->d.normalized(), &newDirection, sampler->next2D());
                //No modification to the fr value as phaseFr/phasePDF = 1
            }
            //Otherwise proceed like in the usual path_mats, but considering the transmittance of the medium as well
            else {
                if (!hitSurface) {
                    // May have hit an environment map!
                    if (scene->getEnvironmentMap()) {
                        Color3f emitted = scene->getEnvironmentMap()->eval(_ray->d);
                        emittedAndReflected += (throughput * emitted * fr / continueProbability);
                    }

                    delete _ray;
                    return emittedAndReflected;
                }

                if (its.mesh->isEmitter()) {
                    Color3f emitted = its.mesh->getEmitter()->eval(its.toLocal(-_ray->d.normalized()), &its);
                    emittedAndReflected += (throughput * emitted * fr / continueProbability);
                }
                
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
                newPoint = its.p;
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

            bounces += 1;
        }

        delete _ray;
        return emittedAndReflected;
    }

    std::string toString() const {
        return "VolPathMats[]";
    }

};

NORI_REGISTER_CLASS(VolPathMats, "vol_path_mats");
NORI_NAMESPACE_END