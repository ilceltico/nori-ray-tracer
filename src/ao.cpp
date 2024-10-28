#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <pcg32.h>

NORI_NAMESPACE_BEGIN

#define MAX(a, b) (((a)>(b))?(a):(b))


using nori::Warp;

class AmbientOcclusionIntegrator : public Integrator {
public:
    AmbientOcclusionIntegrator(const PropertyList &props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);


        /* Sample a new direction on the hemisphere */
        Point2f sample = sampler->next2D();
        Point3f direction = Warp::squareToCosineHemisphere(sample);
        // Point3f direction = Warp::squareToUniformHemisphere(sample);

        /* Orient the hemisphere according to the surface normal */
        Point3f worldDirection = its.shFrame.toWorld(direction);

        Ray3f shadowRay = Ray3f(its.p, worldDirection);
        if (scene->rayIntersect(shadowRay))
            return Color3f(0.0f);

        /* To be used with Uniform Hemisphere */
        // float costheta = worldDirection.normalized().dot(its.shFrame.n.normalized());
        // float coeff = costheta / M_PI * 2 * M_PI;
        // return Color3f(1,1,1) * coeff;

        return Color3f(1,1,1);
    }

    std::string toString() const {
        return "AmbientOcclusionIntegrator[]";
    }

};

NORI_REGISTER_CLASS(AmbientOcclusionIntegrator, "ao");
NORI_NAMESPACE_END