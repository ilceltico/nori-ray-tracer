/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <nori/texture.h>
#include <nori/object.h>
#include <nori/mesh.h>
#include <nori/medium.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Convenience data structure used to pass multiple
 * parameters to the evaluation and sampling routines in \ref BSDF
 */
struct BSDFQueryRecord {
    /// Incident direction (in the local frame)
    Vector3f wi;

    /// Outgoing direction (in the local frame)
    Vector3f wo;

    /// Relative refractive index in the sampled direction
    float eta;

    /// Measure associated with the sample
    EMeasure measure;

    Sampler* sampler = nullptr;

    // Point2f uv = Point2f(std::numeric_limits<float>::quiet_NaN());
    Intersection* its = nullptr;

    /// Create a new record for sampling the BSDF
    BSDFQueryRecord(const Vector3f &wi)
        : wi(wi), eta(1.f), measure(EUnknownMeasure) { }

    /// Create a new record for sampling the BSDF
    BSDFQueryRecord(const Vector3f &wi, Sampler* _sampler)
        : wi(wi), eta(1.f), measure(EUnknownMeasure) { sampler = _sampler; }

    /// Create a new record for querying the BSDF
    BSDFQueryRecord(const Vector3f &wi,
            const Vector3f &wo, EMeasure measure)
        : wi(wi), wo(wo), eta(1.f), measure(measure) { }

    /// Create a new record for querying the BSDF
    BSDFQueryRecord(const Vector3f &wi,
            const Vector3f &wo, EMeasure measure, Sampler* _sampler)
        : wi(wi), wo(wo), eta(1.f), measure(measure) { sampler = _sampler; }

    /// Create a new record for sampling the BSDF
    BSDFQueryRecord(const Vector3f &wi, Intersection* its)
        : wi(wi), eta(1.f), measure(EUnknownMeasure), its(its) { }

    /// Create a new record for sampling the BSDF
    BSDFQueryRecord(const Vector3f &wi, Sampler* _sampler, Intersection* its)
        : wi(wi), eta(1.f), measure(EUnknownMeasure), its(its) { sampler = _sampler; }

    /// Create a new record for querying the BSDF
    BSDFQueryRecord(const Vector3f &wi,
            const Vector3f &wo, EMeasure measure, Intersection* its)
        : wi(wi), wo(wo), eta(1.f), measure(measure), its(its) { }

    /// Create a new record for querying the BSDF
    BSDFQueryRecord(const Vector3f &wi,
            const Vector3f &wo, EMeasure measure, Sampler* _sampler, Intersection* its)
        : wi(wi), wo(wo), eta(1.f), measure(measure), its(its) { sampler = _sampler; }
};

/**
 * \brief Superclass of all bidirectional scattering distribution functions
 */
class BSDF : public NoriObject {
public:

    ~BSDF() {
        delete m_texture;
    }

    virtual void activate() {};

    /**
     * \brief Sample the BSDF and return the importance weight (i.e. the
     * value of the BSDF * cos(theta_o) divided by the probability density
     * of the sample with respect to solid angles).
     *
     * \param bRec    A BSDF query record
     * \param sample  A uniformly distributed sample on \f$[0,1]^2\f$
     *
     * \return The BSDF value divided by the probability density of the sample
     *         sample. The returned value also includes the cosine
     *         foreshortening factor associated with the outgoing direction,
     *         when this is appropriate. A zero value means that sampling
     *         failed.
     */
    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const = 0;

    /**
     * \brief Evaluate the BSDF for a pair of directions and measure
     * specified in \code bRec
     *
     * \param bRec
     *     A record with detailed information on the BSDF query
     * \return
     *     The BSDF value, evaluated for each color channel
     */
    virtual Color3f eval(const BSDFQueryRecord &bRec) const = 0;

    /**
     * \brief Compute the probability of sampling \c bRec.wo
     * (conditioned on \c bRec.wi).
     *
     * This method provides access to the probability density that
     * is realized by the \ref sample() method.
     *
     * \param bRec
     *     A record with detailed information on the BSDF query
     *
     * \return
     *     A probability/density value expressed with respect
     *     to the specified measure
     */

    virtual float pdf(const BSDFQueryRecord &bRec) const = 0;

    /**
     * \brief Return the type of object (i.e. Mesh/BSDF/etc.)
     * provided by this instance
     * */
    EClassType getClassType() const { return EBSDF; }

    /**
     * \brief Return whether or not this BRDF is diffuse. This
     * is primarily used by photon mapping to decide whether
     * or not to store photons on a surface
     */
    virtual bool isDiffuse() const { return false; }

    virtual bool requiresDifferentials() const { return false; }

    /**
     * \brief Returns true only for BRDFs that do not interact with light at all
    */
    virtual bool isNull() const {return false; }



    void addChild(NoriObject *obj) {
        switch (obj->getClassType()) {
            case ETexture:
                if (m_texture)
                    throw NoriException(
                        "BSDF: tried to register multiple Texture instances!");
                m_texture = static_cast<Texture *>(obj);
                break;

            default:
                throw NoriException("BSDF::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    } //#TODO this does not fail at parsing-time in case some BSDFs do not support a texture

protected:
    Texture *m_texture = nullptr;
};

NORI_NAMESPACE_END
