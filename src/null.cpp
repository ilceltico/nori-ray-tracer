
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/texture.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Completely intangible BRDF
 */
class Null : public BSDF {
public:
    Null(const PropertyList &propList) {
    }

    void activate() {
    }

    /// Evaluate the BRDF model
    Color3f eval(const BSDFQueryRecord &bRec) const {
        return Color3f(0.0f);
    }

    /// Compute the density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        return 0.0f;
    }

    /// Draw a a sample from the BRDF model
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        bRec.wo = -bRec.wi;
        return Color3f(1.0f);
    }

    bool isDiffuse() const {
        return false;
    }

    bool isNull() const {
        return true;
    }

    /// Return a human-readable summary
    std::string toString() const {
        return tfm::format(
            "Null[]");
    }

    EClassType getClassType() const { return EBSDF; }

};

NORI_REGISTER_CLASS(Null, "null");
NORI_NAMESPACE_END
