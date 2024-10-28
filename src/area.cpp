#include <nori/emitter.h>
#include <nori/scene.h>
#include <nori/texture.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
public:
    AreaLight(const PropertyList &props) {
        this->m_constant_radiance = props.getColor("radiance", Color3f(100.0f));
    }

    void preprocess(const Scene &scene) {
        if (!m_texture) {
            m_texture = new ConstantTexture(m_constant_radiance);
        }
    }

    Color3f eval(const Vector3f &wo, Intersection *its = nullptr) const {
        if (Frame::cosTheta(wo) < 0) {
                return Color3f(0,0,0);
        }
        BSDFQueryRecord bRec = BSDFQueryRecord(wo, its);
        return m_texture->eval(bRec);
    }

    Color3f sample(const Point3f &x, const Point2f &sample, Point3f &y, Normal3f &yNormal, float &rayEpsilon) const {
        float u, v;
        float pdf = mesh->uniformSquareToUniformSurface(sample, y, yNormal, u , v);
        rayEpsilon = Epsilon;
        Intersection its;
        its.uv = Point2f(u,v);
        BSDFQueryRecord bRec = BSDFQueryRecord(y-x, &its);
        return m_texture->eval(bRec) / pdf;
    }

    float pdf(const Point3f &x, const Point3f &y) const {
        return mesh->getInverseTotalSurfaceArea();
    }

    bool isAssociatedToMesh() const {
        return true;
    }

    void setAssociatedMesh(Mesh *mesh) {
        this->mesh = mesh;
    }

    std::string toString() const {
        return "AreaLight[]";
    }

private:
    Color3f m_constant_radiance;
    Mesh *mesh;
};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END