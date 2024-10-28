
#include <nori/medium.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class Homogeneous : public Medium {
public:
    Homogeneous(const PropertyList &propList) : Medium(propList) {
        m_absorption = propList.getColor("absorption", Color3f(0.5f));
        m_scattering = propList.getColor("scattering", Color3f(0.5f));

        m_extinction = m_absorption + m_scattering;
    }

    Color3f transmittance(const Ray3f &ray, Sampler &sampler) const {
        return Color3f(
            std::exp(-m_extinction[0] * std::min(ray.maxt * ray.d.norm(), std::numeric_limits<float>::max())),
            std::exp(-m_extinction[1] * std::min(ray.maxt * ray.d.norm(), std::numeric_limits<float>::max())),
            std::exp(-m_extinction[2] * std::min(ray.maxt * ray.d.norm(), std::numeric_limits<float>::max()))
            );
    }

    Color3f sample(const Ray3f &ray, Sampler &sampler, Point3f &interactionPoint, bool &sampledMedium) const {
        // Uniformly sample a channel
        int channel = std::min((int)(sampler.next1D() * 3), 2);

        // Sample a distance along the ray
        float dist = -std::log(1 - sampler.next1D()) / m_extinction[channel];
        float t = std::min(dist / ray.d.norm(), ray.maxt);

        // Check if we exceeded the max ray length (i.e. hit a surface)
        sampledMedium = t < ray.maxt;
        if (sampledMedium) {
            interactionPoint = ray(t);
        }

        // Compute the transmittance and sampling density
        Color3f transmittance = Color3f(
                std::exp(-m_extinction[0] * std::min(t, std::numeric_limits<float>::max()) * ray.d.norm()),
                std::exp(-m_extinction[1] * std::min(t, std::numeric_limits<float>::max()) * ray.d.norm()),
                std::exp(-m_extinction[2] * std::min(t, std::numeric_limits<float>::max()) * ray.d.norm())
            );

        // Return weighting factor for scattering from homogeneous medium
        Color3f density = sampledMedium ? (m_extinction * transmittance) : transmittance;
        float pdf = 0;
        for (int i = 0; i < 3; ++i) pdf += density[i];
        pdf *= 1 / 3.0f;
        if (pdf == 0) {
            if (!transmittance.isZero()) {
                throw NoriException("Homogeneous Medium: transmittance is not zero");
            }
            cout << "pdf was zero\n";
            pdf = 1;
        }
        
        if (sampledMedium)
            transmittance *= m_scattering;

        return transmittance / pdf;
    }

    /// Return a human-readable summary
    std::string toString() const {
        return tfm::format(
            "Homogeneous[\n"
            "  absorption = %s\n"
            "  scattering = %s\n"
            "  extinction = %s\n"
            "  phase function = %s\n"
            "]", m_absorption.toString(), m_scattering.toString(), m_extinction.toString(), m_phase ? m_phase->toString() : std::string("null"));
    }

    EClassType getClassType() const { return EMedium; }

private:
    Color3f m_absorption;
    Color3f m_scattering;
    Color3f m_extinction;
};

NORI_REGISTER_CLASS(Homogeneous, "homogeneous");
NORI_NAMESPACE_END
