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

        m_eta = m_intIOR/m_extIOR;
        m_invEta = 1.0f/m_eta;

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
        if (bRec.measure != ESolidAngle || Frame::cosTheta(bRec.wi) == 0)
            return Color3f(0.0f);

        /* Determine the type of interaction */
        bool reflect = Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) > 0;

        // float alpha = (1.2f - 0.2f * sqrt(abs(Frame::cosTheta(bRec.wi)))) * this->m_alpha;
        float alpha = this->m_alpha;

        Vector3f H;
        if (reflect) {

            /* Calculate the reflection half-vector */
            H = (bRec.wo+bRec.wi).normalized();
        } else {
            /* Calculate the transmission half-vector */
            float eta = Frame::cosTheta(bRec.wi) > 0
                ? m_eta : m_invEta;

            H = (bRec.wi + bRec.wo*eta).normalized();
        }

        /* Ensure that the half-vector points into the
           same hemisphere as the macrosurface normal */
        H *= sign(Frame::cosTheta(H));

        /* Evaluate the microfacet normal distribution */
        const float D = Warp::squareToBeckmannPdf(H, alpha) / Frame::cosTheta(H);
        if (D == 0){
            return Color3f(0.0f);
        }

        /* Fresnel factor */
        const float F = fresnel(bRec.wi.dot(H), m_extIOR, m_intIOR);

        /* Smith's shadow-masking function */
        const float G = G1(bRec.wi, H, alpha) * G1(bRec.wo, H, alpha);

        if (reflect) {
            /* Calculate the total amount of reflection */
            float value = F * D * G /
                (4.0f * std::abs(Frame::cosTheta(bRec.wi)) * std::abs(Frame::cosTheta(bRec.wo))); //This last cosine division is needed because Nori's path tracers compute cosines outsides of the BRDFs, and they use the ITS reference frame instead of the microfacet one. This term cancels those cosines.

            return value;
        } else {
            float eta = Frame::cosTheta(bRec.wi) > 0.0f ? m_eta : m_invEta;

            /* Calculate the total amount of transmission */
            float sqrtDenom = bRec.wi.dot(H) + eta * bRec.wo.dot(H);
            float value = ((1 - F) * D * G * eta * eta
                * bRec.wi.dot(H) * bRec.wo.dot(H)) /
                (Frame::cosTheta(bRec.wi) * sqrtDenom * sqrtDenom) / std::abs(Frame::cosTheta(bRec.wo)); //This last cosine division is needed because Nori's path tracers compute cosines outsides of the BRDFs, and they use the ITS reference frame instead of the microfacet one. This term cancels those cosines.

            /* Missing term in the original paper: account for the solid angle
               compression when tracing radiance -- this is necessary for
               bidirectional methods */
            float factor = Frame::cosTheta(bRec.wi) > 0 ? m_invEta : m_eta;
            // float factor = 1.0f;

            return std::abs(value * factor * factor);
        }

        // Color3f result;

        // // Refraction (transmission)
        // if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) < 0) {
        //     Vector3f wi = bRec.wi;
        //     Vector3f wo = bRec.wo;
        //     float eta_i = m_extIOR, eta_o = m_intIOR;
        //     if (Frame::cosTheta(bRec.wi) < 0.0f) {
        //         eta_i = m_intIOR;
        //         eta_o = m_extIOR;
        //     } 
        //     float eta = eta_i / eta_o;
        //     // Vector3f wh = (-eta_i * wi - eta_o * wo).normalized();
        //     Vector3f wh = (wi + wo * eta).normalized();
        //     wh *= sign(Frame::cosTheta(wh));

        //     float D = Warp::squareToBeckmannPdf(wh, this->m_alpha) / Frame::cosTheta(wh);

        //     float F_t = 1 - fresnel(wh.dot(wi), this->m_extIOR, this->m_intIOR);

        //     // float G = G1(wi, wh) * G1(wo, wh); //b?
        //     float G = G1(wi, wh) * G1(wo, wh);

        //     float sqrtDenom = wi.dot(wh) + eta * wo.dot(wh);
        //     float value = (F_t * D * G * eta * eta * wi.dot(wh) * wo.dot(wh)) / (Frame::cosTheta(wi) * sqrtDenom * sqrtDenom);
        //     float factor = 1.0f / eta;
        //     result = abs(value * factor * factor);

        //     // Color3f microfacet = D * F_t * G * eta_o * eta_o;
        //     // microfacet *= wh.dot(wi) * wh.dot(wo);
        //     // microfacet /= (Frame::cosTheta(wi));// * Frame::cosTheta(wo));
        //     // microfacet /= pow(eta_i * wh.dot(wi) + eta_o * wh.dot(wo), 2);
        //     // microfacet = abs(microfacet);

        //     // result = microfacet;

        //     if (!result.isValid()) {
        //         cout << result << "\n";
        //         return Color3f(0.0f);
        //     }
        // }

        // // Reflection
        // else {
        //     Vector3f wi = bRec.wi;
        //     Vector3f wo = bRec.wo;
        //     Vector3f wh = (wi + wo).normalized();
        //     wh *= sign(Frame::cosTheta(wh));

        //     // float D = Warp::squareToBeckmannPdf(wh * sign(Frame::cosTheta(wh)), this->m_alpha);
        //     float D = Warp::squareToBeckmannPdf(wh, this->m_alpha) / Frame::cosTheta(wh);

        //     float F_r = fresnel(wh.dot(wi), this->m_extIOR, this->m_intIOR);

        //     // float G = G1(wi, wh) * G1(wo, wh); //b?
        //     float G = G1(wi, wh) * G1(wo, wh);

        //     Color3f microfacet = D * F_r * G;
        //     microfacet /= (4.0f * abs(Frame::cosTheta(wi)));// * abs(Frame::cosTheta(wo)));// * abs(Frame::cosTheta(wh))); //?

        //     result = microfacet;
        // }
    	
        

        // if (!result.isValid()) {
        //     if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
        //         cout << "Refraction\n";

        //     cout << result.toString() << "\n";
        //     return Color3f(0.0f);
        // }
        // return result;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle) {
            return 0.0f;
        }
        
        Vector3f wh;
        float jacob;

        bool reflect = Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) > 0.0f;

        if (reflect) {
            wh = (bRec.wo + bRec.wi).normalized();
            jacob = 1.0f / (4.0f * bRec.wo.dot(wh));
        }
        else {
            float eta = Frame::cosTheta(bRec.wi) > 0 ? m_eta : m_invEta;

            wh = (bRec.wi + bRec.wo*eta).normalized();

            /* Jacobian of the half-direction mapping */
            float sqrtDenom = bRec.wi.dot(wh) + eta * bRec.wo.dot(wh);
            jacob = (eta*eta * bRec.wo.dot(wh)) / (sqrtDenom*sqrtDenom);
        }

        wh *= sign(Frame::cosTheta(wh));

        float alpha = (1.2f - 0.2f * sqrt(abs(Frame::cosTheta(bRec.wi)))) * this->m_alpha;
        // float alpha = m_alpha;

        float prob = Warp::squareToBeckmannPdf(wh, alpha);

        float Fr = fresnel(bRec.wi.dot(wh), m_extIOR, m_intIOR);

        prob *= reflect ? Fr : (1.0f-Fr);
        return abs(prob * jacob);

    	// Vector3f wi = bRec.wi;
        // Vector3f wo = bRec.wo;
        // // Vector3f wh = (wi + wo).normalized();
        // Vector3f wh;
        // float result;

        // float alpha = (1.2f - 0.2f * sqrt(abs(Frame::cosTheta(bRec.wi)))) * this->m_alpha;

        // // Refraction (transmission)
        // if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0) {
        //     float eta_i = m_extIOR, eta_o = m_intIOR;
        //     if (Frame::cosTheta(bRec.wi) < 0) {
        //         eta_i = m_intIOR;
        //         eta_o = m_extIOR;
        //     } 
        //     wh = (-eta_i * wi - eta_o * wo).normalized();

        //     float D = Warp::squareToBeckmannPdf(wh * sign(Frame::cosTheta(wh)), alpha);

        //     float F_r = fresnel((wh * sign(Frame::cosTheta(wh))).dot(wi), this->m_extIOR, this->m_intIOR);

        //     result = (1 - F_r) * D * eta_o * eta_o * abs(wo.dot(wh)) / pow((eta_i * (wi.dot(wh)) + eta_o * (wo.dot(wh))),2);
        // }

        // // Reflection
        // else {
        //     wh = (wi + wo).normalized();

        //     float D = Warp::squareToBeckmannPdf(wh * sign(Frame::cosTheta(wh)), alpha);

        //     float F_r = fresnel((wh * sign(Frame::cosTheta(wh))).dot(wi), this->m_extIOR, this->m_intIOR);

        //     result = F_r * D / (4 * wo.dot(wh));
        // }

        // return abs(result);
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {

        Point2f sample = _sample;
        float alpha = (1.2f - 0.2f * sqrt(abs(Frame::cosTheta(bRec.wi)))) * this->m_alpha;
        // float alpha = m_alpha;

        Vector3f m = Warp::squareToBeckmann(sample, alpha);
        float microfacetPDF = Warp::squareToBeckmannPdf(m, alpha);
        if (microfacetPDF == 0) {
            return Color3f(0.0f);
        }

        float cosThetaI = bRec.wi.dot(m);

        float F = fresnel(bRec.wi.dot(m), m_extIOR, m_intIOR);
        float etaI = m_extIOR, etaT = m_intIOR;
        float cosThetaIFresnel = bRec.wi.dot(m);
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
        // Inefficient recomputation of cosThetaT. I could simply get it with a parameter from the Fresnel function, or use the paper's formula for the refraction instead, which doesn't need it and is commented right below the refraction formula.

        if (!bRec.sampler) {
            throw NoriException("[Microfacet Dielectric: sample()] Request is missing a sampler!");
        }

        bool sampleReflection = true;
        float sampledValue = bRec.sampler->next1D();
        // float sampledValue = sample.x();
        if (sampledValue > F) {
            sampleReflection = false;
        }
        // cout << sampledValue << "\n";

        Color3f weight = Color3f(1.0f);

        if (sampleReflection) {
            bRec.wo = (2 * (bRec.wi.dot(m)) * m - bRec.wi);
            bRec.eta = 1.0f;

            if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
                return Color3f(0.0f);

            weight *= 1.0f;
        }

        else {
            if (cosThetaT == 0)
                return Color3f(0.0f);

            // float c = bRec.wi.dot(m);
            // float eta = etaFresnel;
            // float eta = cosThetaT<0.0f ? m_eta : m_invEta;

            float eta = m_eta;
            if (cosThetaT < 0.0f)
                eta = 1.0f / eta;
            bRec.wo = m * (bRec.wi.dot(m) * eta + cosThetaT) - bRec.wi * eta;
            // bRec.wo = (eta * c - sign(Frame::cosTheta(bRec.wi)) * sqrt(1 + eta * (c * c - 1))) * m - eta * bRec.wi;

            // float cosTheta2Sqr = 1 - etaFresnel*etaFresnel * (1-cosThetaI*cosThetaI);
            // Frame mFrame = Frame(m);
            // Vector3f wiMframe = mFrame.toLocal(bRec.wi);
            // bRec.wo = Vector3f(-etaFresnel * wiMframe.x(), -etaFresnel * wiMframe.y(), -sqrt(cosTheta2Sqr));
            // if (cosThetaI < 0) {
            //     // cout << "negative costhetaI\n";
            //     bRec.wo[2] = -bRec.wo.z();
            // }
            // bRec.wo = mFrame.toWorld(bRec.wo);
            // bRec.eta = 1.0f / etaFresnel;

            bRec.eta = cosThetaT<0.0f ? m_eta : m_invEta;
            // bRec.eta = 1.0f / bRec.eta;

            if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
                return Color3f(0.0f);


            float factor = cosThetaT < 0 ? m_invEta : m_eta;
            // float factor = cosThetaT<0.0f? m_invEta : m_eta;

            weight *= (factor * factor);
        }

        float evalM = Warp::squareToBeckmannPdf(m, m_alpha) / Frame::cosTheta(m);
        weight *= abs(evalM * G1(bRec.wi, m, m_alpha) * G1(bRec.wo, m, m_alpha) * bRec.wi.dot(m) / (microfacetPDF * Frame::cosTheta(bRec.wi))); 
        //In mitsuba, these distributions are not evaluated with the alpha-rescaling trick, 
        //however if I do so the results are correct with IOR1.7 and alpha0.7, but very incorrect with IOR1.7 alpha0.0001 (using mitsuba itself as reference)
        //So, if use the scaling trick here as well, I get 0.001 mean deviation with IOR1.7 alpha 0.7, and 0.008 with IOR1.7 alpha0.0001

        // if (!sampleReflection) {
        //     cout << bRec.wo.toString() << "\n";
        // }

        bRec.measure = ESolidAngle;

        // cout << bRec.wo.toString() << "\n";

        return weight;

        // Point2f sample = _sample;

        // // cout << fr << "\n";

        // Vector3f wi = bRec.wi;

        // float alpha = (1.2f - 0.2f * sqrt(abs(Frame::cosTheta(bRec.wi)))) * this->m_alpha;

        // Vector3f wh = Warp::squareToBeckmann(sample, alpha);

        // float fr = fresnel(bRec.wi.dot(wh), this->m_extIOR, this->m_intIOR);
        // float etaI = m_extIOR, etaT = m_intIOR;
        // float cosThetaI = bRec.wi.dot(wh);
        // if (cosThetaI < 0.0f) {
        //     std::swap(etaI, etaT);
        //     cosThetaI = -cosThetaI;
        // }
        // float eta = etaI / etaT;
        // float sinThetaTSqr = eta*eta * (1-cosThetaI*cosThetaI);
        // float cosThetaT = std::sqrt(1.0f - sinThetaTSqr);

        // float factor = 1.0f;
    	
        // // Reflection
        // if (sample.x() < fr) {
        //     sample.x() /= fr;

        //     bRec.wo = (2 * abs(wi.dot(wh)) * wh - wi);

        //     bRec.eta = 1.0f;

        //     if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
        //         return Color3f(0.0f);
        // }

        // // Refraction (Transmission)
        // else {
        //     sample.x() = (sample.x() - fr) / (1 - fr);

        //     if (cosThetaT == 0.0f) {
        //         return Color3f(0.0f);
        //     }
            
        //     float c = wi.dot(wh);
        //     bRec.wo = (eta * c - sign(Frame::cosTheta(wi)) * sqrt(1 + eta * (c * c - 1))) * wh - eta * wi;

        //     bRec.eta = cosThetaT < 0.0f ? m_intIOR / m_extIOR : m_extIOR / m_intIOR;

        //     if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
        //         return Color3f(0.0f);

        //     factor = cosThetaT < 0.0f ? m_extIOR / m_intIOR : m_intIOR / m_extIOR;
        // }

        // bRec.measure = ESolidAngle;

        // // if (Frame::cosTheta(wi) * Frame::cosTheta(bRec.wo) <= 0)
        // //     return Color3f(0.0f); ?

        // // float pdfV = pdf(bRec);
        // // float D = Warp::squareToBeckmannPdf(wh, this->m_alpha); //scale alpha
        // float G = G1(wi, wh) * G1(bRec.wo, wh);
        // Color3f result = eval(bRec) * factor * factor * abs(G * wh.dot(wi) / (Frame::cosTheta(wh) * Frame::cosTheta(wi)));
        // return result;

        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * abs(Frame::cosTheta(bRec.wo)) / pdf(bRec);
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
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_eta;
    float m_invEta;

    /// @brief Utility function to compute the shadowing term of the Microfacet BRDF. Assumes inputs to be normalized
    // float G1(const Vector3f &wv, const Vector3f &wh) const {
    //     // float c = wv.dot(wh) / Frame::cosTheta(wv);
    //     int c = sign(wv.dot(wh)) * sign(Frame::cosTheta(wv));
    //     // float chi = c>0 ? 1.0f : 0.0f;
    //     if (c <= 0)
    //         return 0.0f;
    //     float b = 10.0f;
    //     if (Frame::cosTheta(wv) < 1.0f && Frame::cosTheta(wv) > - 1.0f)
    //         b = 1.0f / this->m_alpha * abs(Frame::cosTheta(wv)) / sqrt(1 - pow(Frame::cosTheta(wv), 2));
    //     float last = b>=1.6f ? 1.0f : (3.535f*b+2.181f*b*b)/(1.0f+2.276f*b+2.577f*b*b);
    //     // return chi * last;
    //     return last;
    // }
    float G1(const Vector3f &wv, const Vector3f &wh, const float alpha) const {
        if (wv.dot(wh) * Frame::cosTheta(wv) <= 0)
            return 0.0f;

        float tanTheta = std::abs(Frame::tanTheta(wv));
        if (tanTheta == 0.0f)
            return 1.0f;

        float a = 1.0f / (alpha * tanTheta);
        if (a >= 1.6f)
            return 1.0f;

        /* Use a fast and accurate (<0.35% rel. error) rational
            approximation to the shadowing-masking function */
        float aSqr = a*a;
        return (3.535f * a + 2.181f * aSqr)
                / (1.0f + 2.276f * a + 2.577f * aSqr);
    }
};

NORI_REGISTER_CLASS(MicrofacetDielectric, "microfacet_dielectric");
NORI_NAMESPACE_END
