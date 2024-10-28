#include <nori/emitter.h>
#include <nori/scene.h>
#include <nori/bbox.h>
#include <filesystem/resolver.h>
#include <nori/warp.h>
#include <math.h>

#define MODULOFLOAT(num, max) (fmod((fmod(num, max) + max), max))

NORI_NAMESPACE_BEGIN

class EnvironmentLight : public Emitter {
public:
    EnvironmentLight(const PropertyList &props) {
        m_filename = props.getString("filename");
        filesystem::path filename = getFileResolver()->resolve(m_filename);
        std::string full_filename = filename.str();
        m_mipmap = Mipmap(full_filename);
        m_bitmap = Bitmap(full_filename);

        this->m_radiancePower = props.getColor("radiance_power", Color3f(1.0f));
        this->m_horizontalRotation = props.getFloat("horizontal_rotation", 0.0f);
        m_horizontalRotation = degToRad(m_horizontalRotation);

        this->ratio = 1.0f;
        this->rows = m_mipmap.getRows();
        this->cols = m_mipmap.getCols();
        this->multiplier = std::max(rows, cols);
        if (rows < cols) {
            ratio = (float) cols / (float) rows;
        } else if (cols < rows) {
            ratio = (float) rows / (float) cols;
        }
    }

    Color3f eval(const Vector3f &wo, Intersection *its = nullptr) const {
        Vector3f d = wo.normalized();

        float xLight = ((std::atan2(d.x(), -d.z()) + this->m_horizontalRotation) * INV_PI) * 0.5f;
        xLight = MODULOFLOAT(xLight, 1.0f);
        float yLight = (std::acos(-d.y())) * INV_PI;

        if (rows < cols) {
            yLight = yLight / ratio;
        } else if (cols < rows) {
            xLight = xLight / ratio;
        }

        Point2f p = Point2f(xLight, yLight);

        int y = rows - (int) (p.y() * multiplier) - 1;
        int x = (int) (p.x() * multiplier);

        // cout << x << ", " << y << "\n";
        // cout << rows << ", " << cols << "\n";

        return m_bitmap(y,x) * m_radiancePower;
    }

    Color3f sample(const Point3f &x, const Point2f &sample, Point3f &y, Normal3f &yNormal, float &rayEpsilon) const {
        Point2f lightSample = m_mipmap.squareToHierarchicalSampleWarping(sample);

        int yBitmap = rows - (int) (lightSample.y() * multiplier) - 1;
        int xBitmap = (int) (lightSample.x() * multiplier);

        Color3f emitted = m_bitmap(yBitmap,xBitmap);

        // cout << xBitmap << ", " << yBitmap << "\n";

        float xLight = lightSample.x() * 2 * M_PI - this->m_horizontalRotation;
        float yLight = lightSample.y() * M_PI;

        if (rows < cols) {
            yLight = yLight * ratio;
        } else if (cols < rows) {
            xLight = xLight * ratio;
        }

        Vector3f d = Vector3f(std::sin(xLight)*std::sin(yLight), -std::cos(yLight), -std::cos(xLight)*std::sin(yLight)).normalized();

        yNormal = -d;
        y = x + 2 * d * m_worldRadius;

        float dist2 = (x-y).squaredNorm();
        rayEpsilon = Epsilon / sqrt(dist2);

        float pdf = m_mipmap.squareToHierarchicalSampleWarpingPdf(lightSample) / std::max(std::abs(std::sin(yLight)), Epsilon) / (2*M_PI*M_PI) / ratio;
        // float pdf = m_mipmap.squareToHierarchicalSampleWarpingPdf(lightSample) / (2*M_PI*M_PI);

        Color3f result = emitted * this->m_radiancePower * dist2 / pdf;

        // cout << emitted.toString() << "\n";
        // cout << result.toString() << "\n\n";

        // if (!result.isValid()) {
        //     cout << yNormal.toString() << "\n";
        //     cout << lightSample.x() << "\n";
        //     cout << lightSample.y() << "\n";
        //     m_mipmap.squareToHierarchicalSampleWarping(sample);
        // }

        return result;
    }

    float pdf(const Point3f &x, const Point3f &y) const {
        Vector3f d = (y-x).normalized();

        float xLight = (std::atan2(d.x(), -d.z()) + this->m_horizontalRotation) * INV_PI * 0.5f;
        xLight = MODULOFLOAT(xLight, 1.0f);
        float yLight = (std::acos(-d.y())) * INV_PI;

        float sinTerm = std::max(std::abs(std::sin(yLight)), Epsilon);

        if (rows < cols) {
            yLight = yLight / ratio;
        } else {
            xLight = xLight / ratio;
        }

        Point2f lightSample = Point2f(xLight, yLight);
        return m_mipmap.squareToHierarchicalSampleWarpingPdf(lightSample) / sinTerm / (2*M_PI*M_PI) / ratio;
    }

    bool isAssociatedToMesh() const {
        return false;
    }

    std::string toString() const {
        return tfm::format(
        "EnvironmentlLight[\n"
        "  filename = %s,\n"
        "  worldRadius = %s,\n"
        "  worldCenter = {\n"
        "  %s  }\n"
        "]",
        indent(m_filename),
        m_worldRadius,
        indent(m_worldCenter.toString()));
    }

    void preprocess(const Scene &scene) {
        BoundingBox3f bb = scene.getBoundingBox();
        m_worldRadius = sqrt(bb.getExtents().squaredNorm())/2;
        m_worldCenter = bb.getCenter();
    }

    bool canBeHit() const {
        return true;
    }

    bool isEnvironmentMap() const {
        return true;
    }

    bool dependsOnDistance() const {
        return false; 
    }

private:
    float m_worldRadius;
    Point3f m_worldCenter;

    std::string m_filename;
    Mipmap m_mipmap;
    Bitmap m_bitmap;

    Color3f m_radiancePower;
    float m_horizontalRotation;

    float ratio;
    int rows, cols, multiplier;
};

NORI_REGISTER_CLASS(EnvironmentLight, "envmap");
NORI_NAMESPACE_END