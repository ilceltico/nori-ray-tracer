#include <nori/emitter.h>
#include <nori/scene.h>
#include <nori/bbox.h>

NORI_NAMESPACE_BEGIN

class DirectionalLight : public Emitter {
public:
    DirectionalLight(const PropertyList &props) {
        this->irradiance = props.getColor("irradiance");
        this->direction = props.getVector("direction").normalized();
    }

    Color3f eval(const Vector3f &wo, Intersection *its = nullptr) const {
        throw NoriException("Requested an eval() of a Directional Light source");
    }

    Color3f sample(const Point3f &x, const Point2f &sample, Point3f &y, Normal3f &yNormal, float &rayEpsilon) const {
        yNormal = this->direction;
        y = x - 2 * this->direction * this->m_worldRadius;

        float dist2 = (x-y).squaredNorm();
        rayEpsilon = Epsilon / sqrt(dist2);

        return this->irradiance * dist2;
    }

    float pdf(const Point3f &x, const Point3f &y) const {
        // Similarly to mirrors, the pdf for a directional emitter is a Dirac's delta. In order for it to work with the MIS integrator, it has to return 1.0f
        return 1.0f;
    }

    bool isAssociatedToMesh() const {
        return false;
    }

    std::string toString() const {
        return tfm::format(
        "DirectionalLight[\n"
        "  irradiance = %s,\n"
        "  direction = %s\n"
        "  worldRadius = %s,\n"
        "  worldCenter = {\n"
        "  %s  }\n"
        "]",
        indent(irradiance.toString()),
        indent(direction.toString()),
        m_worldRadius,
        indent(m_worldCenter.toString()));
    }

    void preprocess(const Scene &scene) {
        BoundingBox3f bb = scene.getBoundingBox();
        m_worldRadius = sqrt(bb.getExtents().squaredNorm())/2;
        m_worldCenter = bb.getCenter();

        // Point3f a = Point3f(2,2,2);
        // Point3f b = Point3f(0,0,0);
        // Vector3f prova = a-b;
        // cout << prova.toString() << "\n";
        // cout << prova.squaredNorm() << "\n";
    }

    bool canBeHit() const {
        return false;
    }

    bool dependsOnDistance() const {
        return false; 
    }

private:
    Color3f irradiance;
    Vector3f direction;
    float m_worldRadius;
    Point3f m_worldCenter;

};

NORI_REGISTER_CLASS(DirectionalLight, "directional");
NORI_NAMESPACE_END