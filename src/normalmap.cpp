
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/texture.h>
#include <nori/mesh.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

class Normalmap : public BSDF {
public:
    Normalmap(const PropertyList &propList) {}

    void activate() {
        if (!m_texture) {
            throw NoriException("Normalmap BSDF: no normalmap provided!");
        }
    }

    Frame perturbFrame(const BSDFQueryRecord &bRec) const {
        Frame result;
        Color3f c = m_texture->eval(bRec);
        Normal3f n = Normal3f(c.x(), c.y(), c.z()); //Normalizing here can help in some cases (i.e. final image) but is not consistent with mitsuba

        for (int i=0; i<3; ++i)
            n[i] = 2 * n[i] - 1;

        Frame frame;
        frame.n = bRec.its->shFrame.n;
        frame.s = (bRec.its->dpdu - frame.n * frame.n.dot(bRec.its->dpdu)).normalized();
        frame.t = frame.n.cross(frame.s);

        result.n = frame.toWorld(n).normalized();
        result.s = (bRec.its->dpdu - result.n * result.n.dot(bRec.its->dpdu)).normalized();
        result.t = result.n.cross(result.s);

        return result;
    }

    Color3f eval(const BSDFQueryRecord &bRec) const {
        Intersection perturbed(*bRec.its);
        perturbed.shFrame = perturbFrame(bRec);

        BSDFQueryRecord perturbedQuery(perturbed.toLocal(bRec.its->toWorld(bRec.wi)), perturbed.toLocal(bRec.its->toWorld(bRec.wo)), bRec.measure, &perturbed);

        if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(perturbedQuery.wo) <= 0)
            return Color3f(0.0f);

        return m_baseBSDF->eval(perturbedQuery) / (bRec.its->shFrame.cosTheta(bRec.wo)) * (Frame::cosTheta(perturbedQuery.wo)); //Cosine substitution, since the frame changed
    }

    float pdf(const BSDFQueryRecord &bRec) const {

        Intersection perturbed(*bRec.its);
        perturbed.shFrame = perturbFrame(bRec);

        BSDFQueryRecord perturbedQuery(perturbed.toLocal(bRec.its->toWorld(bRec.wi)), perturbed.toLocal(bRec.its->toWorld(bRec.wo)), bRec.measure, bRec.sampler, &perturbed);
        if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(perturbedQuery.wo) <= 0)
            return 0.0f;

        return m_baseBSDF->pdf(perturbedQuery);
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        Intersection perturbed(*bRec.its);
        perturbed.shFrame = perturbFrame(bRec);

        BSDFQueryRecord perturbedQuery(perturbed.toLocal(bRec.its->toWorld(bRec.wi)), bRec.sampler, &perturbed);
        
        Color3f result = m_baseBSDF->sample(perturbedQuery, sample);
        if (!result.isZero()) {
            bRec.wo = bRec.its->toLocal(perturbed.toWorld(perturbedQuery.wo));
            bRec.eta = perturbedQuery.eta;
            bRec.measure = perturbedQuery.measure;
            if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(perturbedQuery.wo) <= 0)
                return Color3f(0.0f);
        }

        return result;
    }

    bool isDiffuse() const {
        return m_baseBSDF->isDiffuse();
    }

    /// Return a human-readable summary
    std::string toString() const {
        return tfm::format(
            "Normalmap[\n"
            "  normalmap = %s\n"
            "  BaseBSDF = %s\n"
            "]", m_texture ? m_texture->toString() : std::string("null"), 
                m_baseBSDF ? m_baseBSDF->toString() : std::string("null"));
    }

    void addChild(NoriObject *obj) {
        switch (obj->getClassType()) {
            case ETexture:
                if (m_texture)
                    throw NoriException(
                        "Normalmap BSDF: tried to register multiple Texture instances!");
                m_texture = static_cast<Texture *>(obj);
                break;

            case EBSDF:
                if (m_baseBSDF)
                    throw NoriException(
                        "Normalmap BSDF: tried to register multiple BaseBSDF instances!");
                m_baseBSDF = static_cast<BSDF *>(obj);
                break;

            default:
                throw NoriException("Normalmap BSDF::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    } //#TODO this does not fail at parsing-time in case some BSDFs do not support a texture

    EClassType getClassType() const { return EBSDF; }

    bool requiresDifferentials() const {
         return true; 
    }


private:
    BSDF *m_baseBSDF = nullptr;
};

NORI_REGISTER_CLASS(Normalmap, "normalmap");
NORI_NAMESPACE_END
