#pragma once

#include <nori/object.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Superclass of all phase functions
 */
class PhaseFunction : public NoriObject {
public:

    ~PhaseFunction() {
    }

    /**
     * \brief Sample the phase function
     *
     * \param wi    The incident direction
     * \param wo    Will be filled with the outgoing sampled direction
     * \param sample  A uniformly distributed sample on \f$[0,1]^2\f$
     *
     * \return The PDF of the sampled direction
     */
    virtual float sample(const Vector3f &wi, Vector3f *wo, const Point2f &sample) const = 0;

    /**
     * \brief Compute the probability of sampling \c wo
     * (conditioned on \c wi).
     *
     * \param wi    The incident direction
     * \param wo    The outgoing direction
     *
     * \return
     *     A probability/density value
     */

    virtual float pdf(const Vector3f &wi, const Vector3f &wo) const = 0;

    /**
     * \brief Return the type of object (i.e. Mesh/BSDF/etc.)
     * provided by this instance
     * */
    EClassType getClassType() const { return EPhaseFunction; }


};

NORI_NAMESPACE_END
