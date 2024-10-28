#pragma once

#include <nori/object.h>

NORI_NAMESPACE_BEGIN

struct BSDFQueryRecord;

class Texture : public NoriObject {
public:
    virtual Color3f eval(const BSDFQueryRecord &bRec) const = 0;

    EClassType getClassType() const { return ETexture; }

};

class ConstantTexture : public Texture {
public:
    ConstantTexture(const PropertyList propList);
    ConstantTexture(const Color3f albedo);
    Color3f eval(const BSDFQueryRecord &bRec) const;
    std::string toString() const;

private:
    Color3f m_albedo;
};

NORI_NAMESPACE_END