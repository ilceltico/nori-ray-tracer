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

#include <nori/bsdf.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        float cosTheta1 = Frame::cosTheta(bRec.wi);
        float eta = 1.0f;

        /* Solve Fresnel equations */
        float Fr = fresnel(cosTheta1, m_extIOR, m_intIOR);

        /* Refraction */
        if (sample.x() >= Fr) {
            float eta1 = m_extIOR, eta2 = m_intIOR;
            if (cosTheta1 < 0) {
                eta1 = m_intIOR;
                eta2 = m_extIOR;
            }

            eta = eta1 / eta2;
            float cosTheta2Sqr = 1 - eta*eta * (1-cosTheta1*cosTheta1);

            // Refraction in local coordinates
            bRec.wo = Vector3f(-eta * bRec.wi.x(), -eta * bRec.wi.y(), -sqrt(cosTheta2Sqr));

            if (cosTheta1 < 0)
                bRec.wo[2] = -bRec.wo.z();
        }

        /* Reflection */
        else {
            // Reflection in local coordinates
            bRec.wo = Vector3f(
                -bRec.wi.x(),
                -bRec.wi.y(),
                bRec.wi.z()
            );
        }


        bRec.measure = EDiscrete;

        /* Relative index of refraction */
        bRec.eta = 1.0f / eta;

        return Color3f(1.0f / (eta * eta));
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
