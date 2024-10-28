#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

#define MAX(a, b) (((a)>(b))?(a):(b))

class SimpleIntegrator : public Integrator {
public:
    SimpleIntegrator(const PropertyList &props) {
        this->position = props.getPoint("position");
        this->energy = props.getColor("energy");
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);


        /* Check if the point to be rendered is visible from the light source */
        Ray3f shadowRay = Ray3f(its.p, this->position - its.p);
        if (scene->rayIntersect(shadowRay))
            return Color3f(0.0f);


        float costheta = (this->position - its.p).normalized().dot(its.shFrame.n.normalized());
        float coeff = MAX(0, costheta) / (4 * M_PI * M_PI * (this->position - its.p).squaredNorm());

        return this->energy * coeff;
    }

    std::string toString() const {
        return "SimpleIntegrator[]";
    }

private:
    Point3f position;
    Color3f energy;
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END