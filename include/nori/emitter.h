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

#include <nori/object.h>
#include <nori/mesh.h>
#include <nori/texture.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Superclass of all emitters
 */
class Emitter : public NoriObject {
public:

    /**
     * \brief Return the type of object (i.e. Mesh/Emitter/etc.) 
     * provided by this instance
     * */
    EClassType getClassType() const { return EEmitter; }

    /**
     * \brief Returns true if this emitter can be associated to a mesh, 
     * which for now means that it is an AreaEmitter
    */
    virtual bool isAssociatedToMesh() const {return false; }

    /**
     * \brief Sets the associated mesh to the Emitter.
    */
    virtual void setAssociatedMesh(Mesh *mesh) {return; };

    /**
     * \brief Evaluate the Emitter on a direction.
     *
     * \param wo
     *     The requested direction. It's supposed to be in the emitter's reference frame if it's an area emitter, or the World reference frame otherwise.
     * \return
     *     The radiance value on the specified direction
     */
    virtual Color3f eval(const Vector3f &wo, Intersection *its = nullptr) const = 0;

    /**
     * \brief Sample the Emitter and return the emitted radiance divided by the PDF.
     *
     * \param x             The point whose lighting is being evaluated
     * \param sample        A uniformly distributed sample on \f$[0,1]^2\f$
     * \param[out] y        The sampled point
     * \param[out] yNormal  The normal on the sampled point
     *
     * \return The emitted radiance value divided by the probability density of the sampled point.
     */
    virtual Color3f sample(const Point3f &x, const Point2f &sample, Point3f &y, Normal3f &yNormal, float &rayEpsilon) const = 0;

    /**
     * \brief Compute the probability of sampling \c y
     *
     * This method provides access to the probability density that
     * is realized by the \ref sample() method. The behavior of this function is not defined outside of the Emitter.
     *
     * \param x The point whose lighting is being evaluated
     * \param y The requested point on the light source
     *
     * \return
     *     A probability/density value
     */

    virtual float pdf(const Point3f &x, const Point3f &y) const = 0;


    /**
     * \brief Method to be run before using the emitter and after the Scene has been completely generated.
    */
    virtual void preprocess(const Scene &scene) {return; };

    /**
     * \brief Returns true if this emitter has a non-Dirac's delta pdf function
    */
    virtual bool canBeHit() const {return true; }

    /**
     * \brief Returns true only for environment maps
    */
    virtual bool isEnvironmentMap() const {return false; }

    /**
     * \brief Returns true only for emitters that depend on distance (i.e. not for environment maps or directional lights)
    */
    virtual bool dependsOnDistance() const {return true; }

    void addChild(NoriObject *obj) {
        switch (obj->getClassType()) {
            case ETexture:
                if (m_texture)
                    throw NoriException(
                        "Emitter: tried to register multiple Texture instances!");
                m_texture = static_cast<Texture *>(obj);
                break;

            default:
                throw NoriException("Emitter::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    } //#TODO this does not fail at parsing-time in case some BSDFs do not support a texture

protected:
    Texture *m_texture = nullptr;

};

NORI_NAMESPACE_END
