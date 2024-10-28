
#include <nori/object.h>
#include <nori/bitmap.h>
#include <nori/texture.h>
#include <filesystem/resolver.h>
#include <nori/mesh.h>
#include <nori/bsdf.h>

#define MODULO(num, max) (((num % max) + max) % max)

NORI_NAMESPACE_BEGIN

class BitmapTexture : public Texture {
public:

    BitmapTexture(const PropertyList propList) {
        //#TODO support more file types

        m_filename = propList.getString("filename");
        filesystem::path filename = getFileResolver()->resolve(m_filename);
        std::string full_filename = filename.str();
        m_bitmap = Bitmap(full_filename);


        m_scaleU = propList.getFloat("scale_u", 1.0f);
        m_scaleV = propList.getFloat("scale_v", 1.0f);

        m_scaleChannels = propList.getColor("scale_channels", Color3f(1.0f));
    }

    Color3f eval(const BSDFQueryRecord &bRec) const {
        if (bRec.its == nullptr) {
            throw NoriException("BitmapTexture: invalid UV argument. Check the creation of your BSDFQueryRecord.");
        }

        /* Bilinear interpolation, texture repetition */
        int xSize = m_bitmap.rows();
        int ySize = m_bitmap.cols();
        float x = bRec.its->uv.x() * ySize - 0.5f;
        float y = xSize - (bRec.its->uv.y() * xSize) - 0.5f;

        x = x * m_scaleU;
        y = y * m_scaleV;

        int xPos = std::floor(x);
        int yPos = std::floor(y);

        float dx1 = x - xPos, dx2 = 1.0f - dx1,
              dy1 = y - yPos, dy2 = 1.0f - dy1;

        Color3f retrieved = m_bitmap(MODULO(yPos, ySize), MODULO(xPos, xSize)) * dx2 * dy2
             + m_bitmap(MODULO(yPos, ySize), MODULO(xPos+1, xSize)) * dx2 * dy1
             + m_bitmap(MODULO(yPos+1, ySize), MODULO(xPos, xSize)) * dx1 * dy2
             + m_bitmap(MODULO(yPos+1, ySize), MODULO(xPos+1, xSize)) * dx1 * dy1;


        return retrieved.abs() * m_scaleChannels;
    }

    std::string toString() const {
        return tfm::format(
            "BitmapTexture[\n"
            "  albedo = %s\n"
            "]", m_filename);
    }

private:
    std::string m_filename;
    Bitmap m_bitmap;

    float m_scaleU, m_scaleV;

    Color3f m_scaleChannels;
};

NORI_REGISTER_CLASS(BitmapTexture, "bitmap");
NORI_NAMESPACE_END