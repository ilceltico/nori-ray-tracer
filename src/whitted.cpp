#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

#define MAX(a, b) (((a)>(b))?(a):(b))

class WhittedIntegrator : public Integrator {
public:
    WhittedIntegrator(const PropertyList &props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Point3f x = its.p;
        Color3f result;

        /* Diffuse materials */
        if (its.mesh->getBSDF()->isDiffuse()) {
            Color3f emitted = Color3f(0,0,0);
            Color3f reflected = Color3f(0,0,0);

            if (its.mesh->isEmitter()) {
                emitted = its.mesh->getEmitter()->eval(its.toLocal(-ray.d));
            }

            Point2f sample = sampler->next2D();
            if (scene->getEmitters().empty()) {
                return Color3f(0.0f);
            }

            uint32_t emitterIndex = floor(sampler->next1D() * scene->getEmitters().size());

            Point3f samplePosition;
            Normal3f sampleNormal;

            float epsilon;
            Color3f le = scene->getEmitters()[emitterIndex]->sample(x, sample, samplePosition, sampleNormal, epsilon);
            if (sampleNormal.dot(x - samplePosition) < 0) {
                reflected = Color3f(0,0,0);
            }
            else {
                Ray3f shadowRay = Ray3f(x, (samplePosition - x));
                shadowRay.mint = epsilon;
                shadowRay.maxt = 1 - epsilon;
                if (scene->rayIntersect(shadowRay))
                    reflected = Color3f(0.0f);
                else {
                    BSDFQueryRecord query = BSDFQueryRecord(its.toLocal((samplePosition - x).normalized()), its.toLocal(-ray.d), ESolidAngle, sampler, &its);
            
                    Color3f fr = its.mesh->getBSDF()->eval(query);
                    float g = abs(its.shFrame.cosTheta(its.toLocal((samplePosition - x).normalized())) * sampleNormal.dot((x - samplePosition).normalized())) / (x - samplePosition).squaredNorm();
                    reflected = fr.cwiseProduct(le) * g * scene->getEmitters().size();
                }
            }


            result = emitted + reflected;
        }

        /* Specular and refractive materials */
        else {
            BSDFQueryRecord query = BSDFQueryRecord(its.toLocal(-ray.d), sampler, &its);
            Color3f fr = its.mesh->getBSDF()->sample(query, sampler->next2D());
            if (fr.cwiseEqual(Color3f(0.0f)).all()) {
                return fr;
            }
            
            if (sampler->next1D() < 0.95f) {
                result = fr * this->Li(scene, sampler, Ray3f(x, its.toWorld(query.wo))) / 0.95f;
            }
            else {
                result = Color3f(0.0f);
            }
        }

        return result;
    }

    std::string toString() const {
        return "WhittedIntegrator[]";
    }

};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END