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
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class TexturedMicrofacet : public BSDF {
public:
    TexturedMicrofacet(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_constant_kd = propList.getColor("kd", Color3f(0.5f));
    }

    void activate() {
        if (!m_texture) {
            m_texture = new ConstantTexture(m_constant_kd);
        }
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

    	Vector3f wi = bRec.wi;
        Vector3f wo = bRec.wo;
        Vector3f wh = (wi + wo).normalized();

        float D = Warp::squareToBeckmannPdf(wh, this->m_alpha);

        float F = fresnel(wh.dot(wi), this->m_extIOR, this->m_intIOR);

        float G = G1(wi, wh) * G1(wo, wh);

        Color3f microfacet = D * F * G;
        microfacet /= (4 * Frame::cosTheta(wi) * Frame::cosTheta(wo) * Frame::cosTheta(wh));

        Color3f result;
        Color3f kd = m_texture->eval(bRec);
        float ks = 1 - kd.maxCoeff();
        result = kd * M_1_PI;
        result += ks * microfacet;

        // if (!result.isValid()) {
        //     cout << G1(wo, wh) << "\n";
        // }
        return result;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

    	Vector3f wi = bRec.wi;
        Vector3f wo = bRec.wo;
        Vector3f wh = (wi + wo).normalized();


        Color3f kd = m_texture->eval(bRec);
        float ks = 1 - kd.maxCoeff();

        float result = ks * Warp::squareToBeckmannPdf(wh, this->m_alpha) * 1.0f / (wh.dot(wo) * 4) + (1 - ks) * Frame::cosTheta(wo) * M_1_PI;
        return result;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
        if (Frame::cosTheta(bRec.wi) <= 0)
            return Color3f(0.0f);

        Point2f sample = _sample;


        Color3f kd = m_texture->eval(bRec);
        float ks = 1 - kd.maxCoeff();
    	
        if (sample.x() < ks) {
            // Specular case
            sample.x() /= ks;
            Vector3f wh = Warp::squareToBeckmann(sample, this->m_alpha);
            Frame whFrame = Frame(wh);
            Vector3f wiLocal = whFrame.toLocal(bRec.wi);
            Vector3f woLocal = Vector3f(
                -wiLocal.x(),
                -wiLocal.y(),
                wiLocal.z()
            );
            bRec.wo = whFrame.toWorld(woLocal);
        } else {
            // Diffuse case
            sample.x() = (sample.x() - ks) / (1 - ks);
            bRec.wo = Warp::squareToCosineHemisphere(sample);
        }
        bRec.measure = ESolidAngle;
        bRec.eta = 1.0f;

        if (Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        float pdfV = pdf(bRec);
        Color3f result = eval(bRec) * Frame::cosTheta(bRec.wo) / pdfV;
        return result;

        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "TexturedMicrofacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd texture = %s,\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_texture ? m_texture->toString() : std::string("null")
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    Color3f m_constant_kd;

    /// @brief Utility function to compute the shadowing term of the Microfacet BRDF. Assumes inputs to be normalized
    float G1(const Vector3f &wv, const Vector3f &wh) const {
        // float c = wv.dot(wh) / Frame::cosTheta(wv);
        int c = sign(wv.dot(wh)) * sign(Frame::cosTheta(wv));
        // float chi = c>0 ? 1.0f : 0.0f;
        if (c < 0)
            return 0.0f;
        float b = 10.0f;
        if (Frame::cosTheta(wv) <= 1.0f - Epsilon)
            b = 1.0f / this->m_alpha * Frame::cosTheta(wv) / sqrt(1 - pow(Frame::cosTheta(wv), 2));
        float last = b>=1.6f ? 1.0f : (3.535f*b+2.181f*b*b)/(1.0f+2.276f*b+2.577f*b*b);
        // return chi * last;
        return last;
    }
};

NORI_REGISTER_CLASS(TexturedMicrofacet, "textured_microfacet");
NORI_NAMESPACE_END
