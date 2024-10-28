#include <nori/phasefunction.h>
#include <nori/frame.h>
#include <math.h>

NORI_NAMESPACE_BEGIN


class HenyeyGreenstein : public PhaseFunction {
public:
    HenyeyGreenstein(const PropertyList &propList) {
        m_g = propList.getFloat("g", 0.5f);
    }

    float pdf(const Vector3f &wi, const Vector3f &wo) const {
        float cosTheta = wi.dot(wo);
        float denom = 1 + m_g * m_g + 2 * m_g * cosTheta;
        return INV_PI/4 * (1 - m_g * m_g) / (denom * std::sqrt(denom));
    }

    float sample(const Vector3f &wi, Vector3f *wo, const Point2f &sample) const {
        float cosTheta;
        if (std::abs(m_g) < 1e-3)
            cosTheta = 1 - 2 * sample[0];
        else {
            float sqrTerm = (1 - m_g * m_g) / (1 + m_g - 2 * m_g * sample[0]);
            cosTheta = -(1 + m_g * m_g - sqrTerm * sqrTerm) / (2 * m_g);
        }

        // Compute direction _wo_ for Henyey--Greenstein sample
        float sinTheta = std::sqrt(std::max(0.0f, 1.0f - cosTheta * cosTheta));
        float phi = 2 * M_PI * sample[1];

        // Compute the out direction given by the computed spherical coordinates
        Frame f = Frame(wi);
        *wo = sinTheta * std::cos(phi) * f.s + sinTheta * std::sin(phi) * f.t + cosTheta * wi;
        // *wo = sphericalDirection(std::acos(cosTheta), phi);

        return pdf(wi, *wo);
    }

    /// Return a human-readable summary
    std::string toString() const {
        return tfm::format(
            "HenyeyGreenstein[\n"
            "  g = %s\n"
            "]", m_g);
    }

    EClassType getClassType() const { return EPhaseFunction; }

private:
    float m_g;
};

NORI_REGISTER_CLASS(HenyeyGreenstein, "hg");
NORI_NAMESPACE_END
