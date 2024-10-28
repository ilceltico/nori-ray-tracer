#pragma once

#include <nori/phasefunction.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Superclass of all media
 */
class Medium : public NoriObject {
public:

    virtual ~Medium(){}

    Medium(const PropertyList &propList) {
        isInternal = propList.getBoolean("internal", true);
    }

    virtual Color3f transmittance(const Ray3f &ray, Sampler &sampler) const = 0;

    virtual Color3f sample(const Ray3f &ray, Sampler &sampler, Point3f &interactionPoint, bool &sampledMedium) const = 0;

    PhaseFunction* getPhaseFunction() { return m_phase; }

    bool isInternalMedium() const { return isInternal; }

    void addChild(NoriObject *obj) {
        switch (obj->getClassType()) {
            case EPhaseFunction:
                if (m_phase)
                    throw NoriException(
                        "Medium: tried to register multiple Phase Function instances!");
                m_phase = static_cast<PhaseFunction *>(obj);
                break;

            default:
                throw NoriException("Medium::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    }

    EClassType getClassType() const { return EMedium; }

protected:
    PhaseFunction *m_phase = nullptr;
    bool isInternal;
};

NORI_NAMESPACE_END
