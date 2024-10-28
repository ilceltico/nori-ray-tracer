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

class MicrofacetDielectric : public BSDF {
public:
    MicrofacetDielectric(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle)
            return Color3f(0.0f);

        Color3f result;

        // Refraction (transmission)
        if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0) {
            Vector3f wi = bRec.wi;
            Vector3f wo = bRec.wo;
            float eta_i = m_extIOR, eta_o = m_intIOR;
            if (Frame::cosTheta(bRec.wi) < 0) {
                eta_i = m_intIOR;
                eta_o = m_extIOR;
            } 
            Vector3f wh = (-eta_i * wi - eta_o * wo).normalized();

            float D = Warp::squareToBeckmannPdf(wh * sign(Frame::cosTheta(wh)), this->m_alpha);

            float F_t = 1 - fresnel(wi.dot(wh * sign(Frame::cosTheta(wh))), this->m_extIOR, this->m_intIOR);

            // float G = G1(wi, wh) * G1(wo, wh); //b?
            // float G = G1(wi * sign(Frame::cosTheta(wi)), wh) * G1(wo * sign(Frame::cosTheta(wo)), wh);
            float G = G1(wi, wh * sign(Frame::cosTheta(wh))) * G1(wo, wh * sign(Frame::cosTheta(wh)));

            Color3f microfacet = D * F_t * G * eta_o * eta_o;
            microfacet *= wh.dot(wi) * wh.dot(wo);
            microfacet /= (Frame::cosTheta(wi) * Frame::cosTheta(wo));
            microfacet /= pow(eta_i * wh.dot(wi) + eta_o * wh.dot(wo), 2);

            result = microfacet;

            if (!result.isValid()) {
                cout << G << "\n";
                return Color3f(0.0f);
            }
        }

        // Reflection
        else {
            Vector3f wi = bRec.wi;
            Vector3f wo = bRec.wo;
            Vector3f wh = (wi + wo).normalized();

            float D = Warp::squareToBeckmannPdf(wh * sign(Frame::cosTheta(wh)), this->m_alpha);

            float F_r = fresnel(wi.dot(wh * sign(Frame::cosTheta(wh))), this->m_extIOR, this->m_intIOR);

            // float G = G1(wi, wh) * G1(wo, wh); //b?
            // float G = G1(wi * sign(Frame::cosTheta(wi)), wh) * G1(wo * sign(Frame::cosTheta(wo)), wh);
            float G = G1(wi, wh * sign(Frame::cosTheta(wh))) * G1(wo, wh * sign(Frame::cosTheta(wh)));

            Color3f microfacet = D * F_r * G;
            microfacet /= (4 * Frame::cosTheta(wi) * Frame::cosTheta(wo));// * Frame::cosTheta(wh)); ?

            result = microfacet;
        }
    	
        

        // if (!result.isValid()) {
        //     if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
        //         cout << "Refraction\n";

        //     cout << result.toString() << "\n";
        //     return Color3f(0.0f);
        // }
        return result;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
    	Vector3f wi = bRec.wi;
        Vector3f wo = bRec.wo;
        // Vector3f wh = (wi + wo).normalized();
        Vector3f wh;
        float result;

        // Refraction (transmission)
        if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0) {
            float eta_i = m_extIOR, eta_o = m_intIOR;
            if (Frame::cosTheta(bRec.wi) < 0) {
                eta_i = m_intIOR;
                eta_o = m_extIOR;
            } 
            wh = (-eta_i * wi - eta_o * wo).normalized();

            float D = Warp::squareToBeckmannPdf(wh * sign(Frame::cosTheta(wh)), this->m_alpha);

            result = D * abs(Frame::cosTheta(wh)) * eta_o * eta_o * abs(wo.dot(wh)) / pow((eta_i * (wi.dot(wh)) + eta_o * (wo.dot(wh))),2);
        }

        // Reflection
        else {
            wh = (wi + wo).normalized();

            float D = Warp::squareToBeckmannPdf(wh * sign(Frame::cosTheta(wh)), this->m_alpha);
    
            result = D * abs(Frame::cosTheta(wh)) / (4 * wo.dot(wh));
        }

        float Fr = fresnel(wi.dot(wh * sign(Frame::cosTheta(wh))), m_extIOR, m_intIOR);
        if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0) {
            result *= (1.0f - Fr);
        }
        else {
            result *= Fr;
        }

        return result;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {

        Point2f sample = _sample;

        Vector3f wi = bRec.wi;

        Vector3f wh = Warp::squareToBeckmann(sample, this->m_alpha);
        float cosThetaI = bRec.wi.dot(wh);


        float fr = fresnel(wi.dot(wh), this->m_extIOR, this->m_intIOR);
        // cout << fr << "\n";
        float etaI = m_extIOR, etaT = m_intIOR;
        float cosThetaIFresnel = bRec.wi.dot(wh);
        if (cosThetaIFresnel < 0.0f) {
            std::swap(etaI, etaT);
            cosThetaIFresnel = -cosThetaIFresnel;
        }
        float etaFresnel = etaI / etaT;
        float sinThetaTSqr = etaFresnel*etaFresnel * (1-cosThetaIFresnel*cosThetaIFresnel);
        float cosThetaT;
        if (sinThetaTSqr >= 1.0f)
            cosThetaT = 0.0f; 
        else
            cosThetaT = std::sqrt(1.0f - sinThetaTSqr);
        cosThetaT *= -sign(cosThetaI);

        if (!bRec.sampler) {
            throw NoriException("[Microfacet Dielectric: sample()] Request is missing a sampler!");
        }
    	
        // Reflection
        if (bRec.sampler->next1D() < fr) {

            bRec.wo = (2 * wi.dot(wh) * wh - wi);

            if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
                return Color3f(0.0f);

            bRec.eta = 1.0f;
        }

        // Refraction (Transmission)
        else {

            // float eta_i = m_extIOR, eta_o = m_intIOR;
            // if (Frame::cosTheta(bRec.wi) < 0) {
            //     eta_i = m_intIOR;
            //     eta_o = m_extIOR;
            // } 
            // float eta = eta_i / eta_o;
            // float c = wi.dot(wh);
            // bRec.wo = (eta * c - sign(Frame::cosTheta(wi)) * sqrt(1 + eta * (c * c - 1))) * wh - eta * wi;
            
            float eta = m_intIOR/m_extIOR;
            if (cosThetaI > 0.0f)
                eta = 1.0f / eta;
            bRec.wo = wh * (bRec.wi.dot(wh) * eta + cosThetaT) - bRec.wi * eta;



            bRec.eta = 1.0f / eta;

            if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0) {
                // cout <<  wi.toString() << "wwww\n";
                return Color3f(0.0f);
            }
        }

        bRec.measure = ESolidAngle;

        float pdfV = pdf(bRec);
        Color3f result = eval(bRec) * abs(Frame::cosTheta(bRec.wo)) / pdfV; //?

        // cout << result.toString() << "\n";

        // if (pdfV <= 0) {
        //     cout << pdfV << "\n";
        //     return Color3f(0.0f);
        // }

        if (!result.isValid()) {
            cout << pdfV << "\n";
            return Color3f(0.0f);
        }


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
            "MicrofacetDielectric[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;

    /// @brief Utility function to compute the shadowing term of the Microfacet BRDF. Assumes inputs to be normalized
    float G1(const Vector3f &wv, const Vector3f &wh) const {
        // float c = wv.dot(wh) / Frame::cosTheta(wv);
        int c = sign(wv.dot(wh)) * sign(Frame::cosTheta(wv));
        // float chi = c>0 ? 1.0f : 0.0f;
        if (c < 0)
            return 0.0f;
        float b = 10.0f;
        if (Frame::cosTheta(wv) <= 1.0f - Epsilon && Frame::cosTheta(wv) >= - 1.0f + Epsilon)
            b = 1.0f / this->m_alpha * Frame::cosTheta(wv) / sqrt(1 - pow(Frame::cosTheta(wv), 2));
        float last = b>=1.6f ? 1.0f : (3.535f*b+2.181f*b*b)/(1.0f+2.276f*b+2.577f*b*b);
        // return chi * last;
        return last;
    }
};

NORI_REGISTER_CLASS(MicrofacetDielectric, "microfacet_dielectric");
NORI_NAMESPACE_END
