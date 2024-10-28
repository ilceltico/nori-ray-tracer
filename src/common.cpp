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

#include <nori/object.h>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <filesystem/resolver.h>
#include <iomanip>
#include <nori/timer.h>
#include <fstream>

#if defined(PLATFORM_LINUX)
#include <malloc.h>
#endif

#if defined(PLATFORM_WINDOWS)
#include <windows.h>
#endif

#if defined(PLATFORM_MACOS)
#include <sys/sysctl.h>
#endif

NORI_NAMESPACE_BEGIN

std::string indent(const std::string &string, int amount) {
    /* This could probably be done faster (it's not
       really speed-critical though) */
    std::istringstream iss(string);
    std::ostringstream oss;
    std::string spacer(amount, ' ');
    bool firstLine = true;
    for (std::string line; std::getline(iss, line); ) {
        if (!firstLine)
            oss << spacer;
        oss << line;
        if (!iss.eof())
            oss << endl;
        firstLine = false;
    }
    return oss.str();
}

bool endsWith(const std::string &value, const std::string &ending) {
    if (ending.size() > value.size())
        return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

std::string toLower(const std::string &value) {
    std::string result;
    result.resize(value.size());
    std::transform(value.begin(), value.end(), result.begin(), ::tolower);
    return result;
}

bool toBool(const std::string &str) {
    std::string value = toLower(str);
    if (value == "false")
        return false;
    else if (value == "true")
        return true;
    else
        throw NoriException("Could not parse boolean value \"%s\"", str);
}

int toInt(const std::string &str) {
    char *end_ptr = nullptr;
    int result = (int) strtol(str.c_str(), &end_ptr, 10);
    if (*end_ptr != '\0')
        throw NoriException("Could not parse integer value \"%s\"", str);
    return result;
}

unsigned int toUInt(const std::string &str) {
    char *end_ptr = nullptr;
    unsigned int result = (int) strtoul(str.c_str(), &end_ptr, 10);
    if (*end_ptr != '\0')
        throw NoriException("Could not parse integer value \"%s\"", str);
    return result;
}

float toFloat(const std::string &str) {
    char *end_ptr = nullptr;
    float result = (float) strtof(str.c_str(), &end_ptr);
    if (*end_ptr != '\0')
        throw NoriException("Could not parse floating point value \"%s\"", str);
    return result;
}

Eigen::Vector3f toVector3f(const std::string &str) {
    std::vector<std::string> tokens = tokenize(str);
    if (tokens.size() != 3)
        throw NoriException("Expected 3 values");
    Eigen::Vector3f result;
    for (int i=0; i<3; ++i)
        result[i] = toFloat(tokens[i]);
    return result;
}

std::vector<std::string> tokenize(const std::string &string, const std::string &delim, bool includeEmpty) {
    std::string::size_type lastPos = 0, pos = string.find_first_of(delim, lastPos);
    std::vector<std::string> tokens;

    while (lastPos != std::string::npos) {
        if (pos != lastPos || includeEmpty)
            tokens.push_back(string.substr(lastPos, pos - lastPos));
        lastPos = pos;
        if (lastPos != std::string::npos) {
            lastPos += 1;
            pos = string.find_first_of(delim, lastPos);
        }
    }

    return tokens;
}

std::string timeString(double time, bool precise) {
    if (std::isnan(time) || std::isinf(time))
        return "inf";

    std::string suffix = "ms";
    if (time > 1000) {
        time /= 1000; suffix = "s";
        if (time > 60) {
            time /= 60; suffix = "m";
            if (time > 60) {
                time /= 60; suffix = "h";
                if (time > 12) {
                    time /= 12; suffix = "d";
                }
            }
        }
    }

    std::ostringstream os;
    os << std::setprecision(precise ? 4 : 1)
       << std::fixed << time << suffix;

    return os.str();
}

std::string memString(size_t size, bool precise) {
    double value = (double) size;
    const char *suffixes[] = {
        "B", "KiB", "MiB", "GiB", "TiB", "PiB"
    };
    int suffix = 0;
    while (suffix < 5 && value > 1024.0f) {
        value /= 1024.0f; ++suffix;
    }

    std::ostringstream os;
    os << std::setprecision(suffix == 0 ? 0 : (precise ? 4 : 1))
       << std::fixed << value << " " << suffixes[suffix];

    return os.str();
}

filesystem::resolver *getFileResolver() {
    static filesystem::resolver *resolver = new filesystem::resolver();
    return resolver;
}

Color3f Color3f::toSRGB() const {
    Color3f result;

    for (int i=0; i<3; ++i) {
        float value = coeff(i);

        if (value <= 0.0031308f)
            result[i] = 12.92f * value;
        else
            result[i] = (1.0f + 0.055f)
                * std::pow(value, 1.0f/2.4f) -  0.055f;
    }

    return result;
}

Color3f Color3f::toLinearRGB() const {
    Color3f result;

    for (int i=0; i<3; ++i) {
        float value = coeff(i);

        if (value <= 0.04045f)
            result[i] = value * (1.0f / 12.92f);
        else
            result[i] = std::pow((value + 0.055f)
                * (1.0f / 1.055f), 2.4f);
    }

    return result;
}

bool Color3f::isValid() const {
    for (int i=0; i<3; ++i) {
        float value = coeff(i);
        if (value < 0 || !std::isfinite(value))
            return false;
    }
    return true;
}

float Color3f::getLuminance() const {
    return coeff(0) * 0.212671f + coeff(1) * 0.715160f + coeff(2) * 0.072169f;
}

Transform::Transform(const Eigen::Matrix4f &trafo)
    : m_transform(trafo), m_inverse(trafo.inverse()) { }

std::string Transform::toString() const {
    std::ostringstream oss;
    oss << m_transform.format(Eigen::IOFormat(4, 0, ", ", ";\n", "", "", "[", "]"));
    return oss.str();
}

Transform Transform::operator*(const Transform &t) const {
    return Transform(m_transform * t.m_transform,
        t.m_inverse * m_inverse);
}

Vector3f sphericalDirection(float theta, float phi) {
    float sinTheta, cosTheta, sinPhi, cosPhi;

    sincosf(theta, &sinTheta, &cosTheta);
    sincosf(phi, &sinPhi, &cosPhi);

    return Vector3f(
        sinTheta * cosPhi,
        sinTheta * sinPhi,
        cosTheta
    );
}

Point2f sphericalCoordinates(const Vector3f &v) {
    Point2f result(
        std::acos(v.z()),
        std::atan2(v.y(), v.x())
    );
    if (result.y() < 0)
        result.y() += 2*M_PI;
    return result;
}

void coordinateSystem(const Vector3f &a, Vector3f &b, Vector3f &c) {
    if (std::abs(a.x()) > std::abs(a.y())) {
        float invLen = 1.0f / std::sqrt(a.x() * a.x() + a.z() * a.z());
        c = Vector3f(a.z() * invLen, 0.0f, -a.x() * invLen);
    } else {
        float invLen = 1.0f / std::sqrt(a.y() * a.y() + a.z() * a.z());
        c = Vector3f(0.0f, a.z() * invLen, -a.y() * invLen);
    }
    b = c.cross(a);
}

float fresnel(float cosThetaI, float extIOR, float intIOR) {
    float etaI = extIOR, etaT = intIOR;

    if (extIOR == intIOR)
        return 0.0f;

    /* Swap the indices of refraction if the interaction starts
       at the inside of the object */
    if (cosThetaI < 0.0f) {
        std::swap(etaI, etaT);
        cosThetaI = -cosThetaI;
    }

    /* Using Snell's law, calculate the squared sine of the
       angle between the normal and the transmitted ray */
    float eta = etaI / etaT,
          sinThetaTSqr = eta*eta * (1-cosThetaI*cosThetaI);

    if (sinThetaTSqr > 1.0f)
        return 1.0f;  /* Total internal reflection! */

    float cosThetaT = std::sqrt(1.0f - sinThetaTSqr);

    float Rs = (etaI * cosThetaI - etaT * cosThetaT)
             / (etaI * cosThetaI + etaT * cosThetaT);
    float Rp = (etaT * cosThetaI - etaI * cosThetaT)
             / (etaT * cosThetaI + etaI * cosThetaT);

    return (Rs * Rs + Rp * Rp) / 2.0f;
}

// float fresnelConductor(float cosThetaI, float eta, float k) {
//     /* Modified from "Optics" by K.D. Moeller, University Science Books, 1988 */

//     float cosThetaI2 = cosThetaI*cosThetaI,
//           sinThetaI2 = 1-cosThetaI2,
//           sinThetaI4 = sinThetaI2*sinThetaI2;

//     float temp1 = eta*eta - k*k - sinThetaI2,
//           a2pb2 = std::sqrt(std::max(0.0f, temp1*temp1 + 4*k*k*eta*eta)),
//           a     = std::sqrt(std::max(0.0f, 0.5f * (a2pb2 + temp1)));

//     float term1 = a2pb2 + cosThetaI2,
//           term2 = 2*a*cosThetaI;

//     float Rs2 = (term1 - term2) / (term1 + term2);

//     float term3 = a2pb2*cosThetaI2 + sinThetaI4,
//           term4 = term2*sinThetaI2;

//     float Rp2 = Rs2 * (term3 - term4) / (term3 + term4);

//     return 0.5f * (Rp2 + Rs2);
// }

// Based on PBRTv3
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
Color3f fresnelConductor(float cosThetaI, const Color3f &eta, const Color3f &k) {
    cosThetaI = clamp(cosThetaI, -1.0f, 1.0f);

    float cosThetaI2 = cosThetaI * cosThetaI;
    float sinThetaI2 = 1. - cosThetaI2;
    Color3f eta2 = eta * eta;
    Color3f etak2 = k * k;

    Color3f t0 = eta2 - etak2 - sinThetaI2;
    Color3f a2plusb2 = sqrt(t0 * t0 + 4 * eta2 * etak2);
    Color3f t1 = a2plusb2 + cosThetaI2;
    Color3f a = sqrt(0.5f * (a2plusb2 + t0));
    Color3f t2 = (float)2 * cosThetaI * a;
    Color3f Rs = (t1 - t2) / (t1 + t2);

    Color3f t3 = cosThetaI2 * a2plusb2 + sinThetaI2 * sinThetaI2;
    Color3f t4 = t2 * sinThetaI2;
    Color3f Rp = Rs * (t3 - t4) / (t3 + t4);

    return 0.5 * (Rp + Rs);
}

//Based on PBRTv3
template <typename Predicate>
int findInterval(int size, const Predicate &pred) {
    int first = 0, len = size;
    while (len > 0) {
        int half = len >> 1, middle = first + half;
        // Bisect range based on value of _pred_ at _middle_
        if (pred(middle)) {
            first = middle + 1;
            len -= half + 1;
        } else
            len = half;
    }
    return clamp(first - 1, 0, size - 2);
}


//Based on PBRTv3
inline float Lerp(float t, float v1, float v2) { return (1 - t) * v1 + t * v2; }
//Based on PBRTv3
float interpolateConductorSamples(const float *lambda, const float *vals, int n,
                                 float l) {
    for (int i = 0; i < n - 1; ++i) {
        if (lambda[i + 1] < lambda[i]) {
            throw NoriException("interpolateConductorSamples: unsorted samples");
        }
    }
    if (l <= lambda[0]) return vals[0];
    if (l >= lambda[n - 1]) return vals[n - 1];
    int offset = findInterval(n, [&](int index) { return lambda[index] <= l; });
    if (!(l >= lambda[offset] && l <= lambda[offset + 1])) {
        throw NoriException("interpolateConductorSamples: inconsistent values");
    }
    float t = (l - lambda[offset]) / (lambda[offset + 1] - lambda[offset]);
    return Lerp(t, vals[offset], vals[offset + 1]);
}
//Based on PBRTv3
bool spectrumSamplesSorted(const float *lambda, const float *vals, int n) {
    for (int i = 0; i < n - 1; ++i)
        if (lambda[i] > lambda[i + 1]) return false;
    return true;
}
//Based on PBRTv3
void sortSpectrumSamples(float *lambda, float *vals, int n) {
    std::vector<std::pair<float, float>> sortVec;
    sortVec.reserve(n);
    for (int i = 0; i < n; ++i)
        sortVec.push_back(std::make_pair(lambda[i], vals[i]));
    std::sort(sortVec.begin(), sortVec.end());
    for (int i = 0; i < n; ++i) {
        lambda[i] = sortVec[i].first;
        vals[i] = sortVec[i].second;
    }
}

//Based on PBRTv3
Color3f fromSampledConductor(const float *lambda, const float *v, int n) {
    // Sort samples if unordered, use sorted for returned spectrum
    if (!spectrumSamplesSorted(lambda, v, n)) {
        std::vector<float> slambda(&lambda[0], &lambda[n]);
        std::vector<float> sv(&v[0], &v[n]);
        sortSpectrumSamples(&slambda[0], &sv[0], n);
        return fromSampledConductor(&slambda[0], &sv[0], n);
    }
    float xyz[3] = {0, 0, 0};
    for (int i = 0; i < nCIESamples; ++i) {
        float val = interpolateConductorSamples(lambda, v, n, CIE_lambda[i]);
        xyz[0] += val * CIE_X[i];
        xyz[1] += val * CIE_Y[i];
        xyz[2] += val * CIE_Z[i];
    }
    float scale = float(CIE_lambda[nCIESamples - 1] - CIE_lambda[0]) /
                    float(CIE_Y_integral * nCIESamples);
    xyz[0] *= scale;
    xyz[1] *= scale;
    xyz[2] *= scale;

    Color3f result;
    result[0] = 3.240479f * xyz[0] - 1.537150f * xyz[1] - 0.498535f * xyz[2];
    result[1] = -0.969256f * xyz[0] + 1.875991f * xyz[1] + 0.041556f * xyz[2];
    result[2] = 0.055648f * xyz[0] - 0.204043f * xyz[1] + 1.057311f * xyz[2];

    return result;
}

void readSampledConductor(std::string filename, std::vector<float> &lambdas, std::vector<float> &values) {
    std::ifstream infile(filename);
    
    if (infile.bad() || infile.fail())
        throw NoriException(tfm::format("ReadSampledConductor: impossible to read the file ", filename).c_str());

    std::string line;
    while (true) {
        if (!std::getline(infile, line))
            break;
        std::string::size_type start = line.find_first_not_of(" \t\r\n");
        std::string::size_type end = line.find_last_not_of(" \t\r\n");
        line = line.substr(start == std::string::npos ? 0 : start, end == std::string::npos ? line.length() - 1 : end - start + 1); //Trim
        if (line.length() == 0 || line[0] == '#') //Skip empty lines and comments
            continue;
        std::istringstream iss(line);
        float lambda, value;
        if (!(iss >> lambda >> value))
            break;
        lambdas.push_back(lambda);
        values.push_back(value);
    }

    if (lambdas.size() == 0)
        throw NoriException(tfm::format("ReadSampledConductor: couldn't find any value in the file ", filename).c_str());
}


void printProgress(double current, double max, Timer &timer) {
    double percentage = current/max;
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    double etaMillis = timer.elapsed() / percentage * (1 - percentage);
    printf("\rElapsed: %s, ETA: %s [%.*s%*s] %3d%%, %d/%d %s",  timer.elapsedString().c_str(), timeString(etaMillis, false).c_str(), lpad, PBSTR, rpad, "", val, (int) current, (int) max, "       ");
    fflush(stdout);
}

NORI_NAMESPACE_END
