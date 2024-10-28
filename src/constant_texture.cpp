
#include <nori/object.h>
#include <nori/bitmap.h>
#include <nori/texture.h>
#include <nori/mesh.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

ConstantTexture::ConstantTexture(const PropertyList propList) {
    m_albedo = propList.getColor("albedo", Color3f(0.5f));
}

ConstantTexture::ConstantTexture(const Color3f albedo) {
    m_albedo = albedo;
}

Color3f ConstantTexture::eval(const BSDFQueryRecord &bRec) const {
    return m_albedo;
}

std::string ConstantTexture::toString() const {
    return tfm::format(
        "ConstantTexture[\n"
        "  albedo = %s\n"
        "]", m_albedo.toString());
}

NORI_REGISTER_CLASS(ConstantTexture, "constant");
NORI_NAMESPACE_END