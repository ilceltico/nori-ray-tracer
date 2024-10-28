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

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>
#include <nori/bitmap.h>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) {
    Point2f result = Point2f();
    
    if (sample.x() >= 0.5f)
        result(0) = 1.0f - sqrt(2 * (1.0f - sample.x()));
    else
        result(0) = sqrt(2 * sample.x()) - 1.0f;

    if (sample.y() >= 0.5f)
        result(1) = 1.0f - sqrt(2 * (1.0f - sample.y()));
    else
        result(1) = sqrt(2 * sample.y()) - 1.0f;

    return result;
}

float Warp::squareToTentPdf(const Point2f &p) {
    if ((p.array() < -1.0f).any() || (p.array() > 1.0f).any())
        return 0.0f;
    
    return (1.0f - abs(p.x())) * (1.0f - abs(p.y()));
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    Point2f result = Point2f();
    float r = sqrt(sample.x());
    float phi = 2 * M_PI * sample.y();

    result(0) = r * cos(phi);
    result(1) = r * sin(phi);

    return result;
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    if ((p.x() * p.x() + p.y() * p.y()) > 1.0f)
        return 0.0f;

    return 1 / M_PI;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    Vector3f result = Vector3f();
    float phi = 2 * M_PI * sample.x();
    float costheta = 1 - 2 * sample.y();
    
    float sintheta = sqrt(1 - costheta * costheta);

    result(0) = cos(phi) * sintheta;
    result(1) = sin(phi) * sintheta;
    result(2) = costheta;

    return result;
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    if (abs(v.dot(v) - 1.0f) > Epsilon)
        return 0.0f;

    return 1 / (4 * M_PI);
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    Vector3f result = Vector3f();
    float phi = 2 * M_PI * sample.x();
    float costheta = 1 - sample.y();
    
    float sintheta = sqrt(1 - costheta * costheta);

    result(0) = cos(phi) * sintheta;
    result(1) = sin(phi) * sintheta;
    result(2) = costheta;

    return result;
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    if (abs(v.dot(v) - 1.0f) > Epsilon)
        return 0.0f;

    if (v.z() < 0)
        return 0.0f;

    return 1 / (2 * M_PI);
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    Vector3f result = Vector3f();
    float phi = 2 * M_PI * sample.x();
    float sintheta = sqrt(sample.y());
    float costheta = sqrt(1 - sample.y());    
    

    result(0) = cos(phi) * sintheta;
    result(1) = sin(phi) * sintheta;
    result(2) = costheta;

    return result;
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    if (abs(v.dot(v) - 1.0f) > Epsilon)
        return 0.0f;

    if (v.z() < 0)
        return 0.0f;

    return v.z() / M_PI;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    Vector3f result = Vector3f();
    float phi = 2 * M_PI * sample.x();
    float tantheta = sqrt(- alpha * alpha * log(1 - sample.y()));
    
    float costheta2 = 1.0f / (1 + tantheta * tantheta);
    float sintheta = sqrt(1 - costheta2);

    result(0) = cos(phi) * sintheta;
    result(1) = sin(phi) * sintheta;
    result(2) = sqrt(costheta2);

    return result;
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    if (abs(m.dot(m) - 1.0f) > Epsilon)
        return 0.0f;

    if (m.z() < 0)
        return 0.0f;

    float tan2 = (m.x() * m.x() + m.y() * m.y()) / (m.z() * m.z());
    if (tan2 == INFINITY)
        return 0.0f;
    float cos3 = pow(m.z(), 3);
    float alpha2 = alpha * alpha;

    float result = pow(M_E, -tan2 / alpha2) / (M_PI * alpha2 * cos3);

    if (result < 1e-20f) //Numerical issues prevention, from mitsuba
        result = 0;

    return result;
}

Point2f Warp::squareToHierarchicalSampleWarping(const Point2f &sample, const Mipmap &mipmap) {
    return mipmap.squareToHierarchicalSampleWarping(sample);
}

float Warp::squareToHierarchicalSampleWarpingPdf(const Point2f &p, const Mipmap &mipmap) {
    return mipmap.squareToHierarchicalSampleWarpingPdf(p);
}


void Mipmap::buildMipmap(std::string filename) {
    Bitmap bitmap(filename);
    int rows = bitmap.rows();
    int cols = bitmap.cols();

    if (!((rows & (rows - 1)) == 0))
        throw NoriException("X axis is not a power of 2");

    if (!((cols & (cols - 1)) == 0))
        throw NoriException("Y axis is not a power of 2");

    // if (rows < 2)
    //     throw NoriException("Specified image is too small");

    MatrixXf luminance(rows, cols);
    Color3f color;
    for (int i=0; i<rows; i++) {
        for (int j=0; j<cols; j++){
            color = bitmap(i,j);
            // EXR files are linear RGB (not gamma-encoded), thus using standard coefficients
            luminance(i,j) = 0.2126 * color(0) + 0.7152 * color(1) + 0.0722 * color(2);
        }
    }

    luminance = luminance / luminance.sum();

    mipmap = std::vector<MatrixXf>();
    mipmap.push_back(luminance);
    
    int currentDimRows = rows > 1 ? rows/2 : rows;
    int currentDimCols = cols > 1 ? cols/2 : cols;
    int currentDim = std::min(currentDimRows, currentDimCols);
    while (currentDim > (rows==cols ? 1 : 0)) {
        MatrixXf luminanceReduced(currentDimRows, currentDimCols);
        for (int i=0; i<currentDimRows; i++){
            for (int j=0; j<currentDimCols; j++){
                luminanceReduced(i,j) = mipmap.back()(i*2,j*2) + mipmap.back()(i*2+1,j*2) + mipmap.back()(i*2,j*2+1) + mipmap.back()(i*2+1,j*2+1);            
            }
        }

        mipmap.push_back(luminanceReduced);
        currentDimRows = currentDimRows/2;
        currentDimCols = currentDimCols/2;
        currentDim /= 2;
    }

    currentDim = std::max(currentDimCols, currentDimRows);

    while (currentDim > 1) {
        MatrixXf luminanceReduced(currentDimRows, currentDimCols);
        for (int i=0; i<currentDimRows; i++){
            for (int j=0; j<currentDimCols; j++){
                if (currentDimRows > 1)
                    luminanceReduced(i,j) = mipmap.back()(i*2,j*2) + mipmap.back()(i*2+1,j*2);    
                else        
                    luminanceReduced(i,j) = mipmap.back()(i*2,j*2) + mipmap.back()(i*2,j*2+1);    
            }
        }

        mipmap.push_back(luminanceReduced);
        currentDimRows = currentDimRows/2;
        currentDimCols = currentDimCols/2;
        currentDim /= 2;
    }
}


Point2f Mipmap::squareToHierarchicalSampleWarping(const Point2f &sample) const {
    float sampleX = sample.x();
    float sampleY = sample.y();
    
    int xDisplacement = 0;
    int yDisplacement = 0;
     
    for (int i=mipmap.size()-1; i>=0; i--){
        int cols = mipmap[i].cols();
        int rows = mipmap[i].rows();

        xDisplacement = rows>1 ? 2*xDisplacement : xDisplacement;
        yDisplacement = cols>1 ? 2*yDisplacement : yDisplacement;
        
        float topLuminance = cols>1 ? mipmap[i](yDisplacement,xDisplacement) + mipmap[i](yDisplacement,xDisplacement+1) : mipmap[i](yDisplacement,xDisplacement);
        float bottomLuminance = rows>1 ? mipmap[i](yDisplacement+1,xDisplacement) + mipmap[i](yDisplacement+1,xDisplacement+1) : 0.0f;
        float topNormalizedLuminance = (topLuminance) / (topLuminance + bottomLuminance);
        float bottomNormalizedLuminance = 1 - topNormalizedLuminance;
        float topLeftNormalizedLuminance = cols>1 ? mipmap[i](yDisplacement,xDisplacement) / (topLuminance) : 1.0f;
        float bottomLeftNormalizedLuminance = rows>1 ? (cols>1 ? mipmap[i](yDisplacement+1,xDisplacement) / (bottomLuminance) : 1.0f) : 0.0f;

        if (rows>1 && sampleY < bottomNormalizedLuminance) {
            sampleY = sampleY / bottomNormalizedLuminance;
            yDisplacement += 1;
            if (sampleX < bottomLeftNormalizedLuminance)
                sampleX = sampleX / bottomLeftNormalizedLuminance;
            else {
                sampleX = (sampleX-bottomLeftNormalizedLuminance) / (1-bottomLeftNormalizedLuminance);
                xDisplacement += 1;
            }
        } else {
            sampleY = (sampleY-bottomNormalizedLuminance) / topNormalizedLuminance;
            if (sampleX < topLeftNormalizedLuminance)
                sampleX = sampleX / topLeftNormalizedLuminance;
            else {
                sampleX = (sampleX-topLeftNormalizedLuminance) / (1-topLeftNormalizedLuminance);
                xDisplacement += 1;
            }
        }

    }

    sampleY = 1.0f - sampleY; //Makes sampleY \in (0,1], which is needed for the inverse mapping in the result.
    
    int divider = std::max(mipmap[0].rows(), mipmap[0].cols());

    Point2f result = Point2f(
                ((float) xDisplacement + sampleX) / (float) divider,
                ((float) mipmap[0].rows() - (float) yDisplacement - sampleY) / divider
                );

    return result;
}


float Mipmap::squareToHierarchicalSampleWarpingPdf(const Point2f &p) const {
    if ((p.array() < 0).any() || (p.array() >= 1).any())
        return 0.0f;

    if (mipmap[0].rows() != mipmap[0].cols()) {
        int indexP = mipmap[0].rows() < mipmap[0].cols() ? 1 : 0;
        float ratio = mipmap[0].rows() < mipmap[0].cols() ? (float) mipmap[0].rows() / (float) mipmap[0].cols() : (float) mipmap[0].cols() / (float) mipmap[0].rows();
        if (p[indexP] >= ratio){
            return 0.0f;
        }
    }

    int multiplier = std::max(mipmap[0].rows(), mipmap[0].cols());

    int y = mipmap[0].rows() - (int) (p.y() * multiplier) - 1;
    int x = (int) (p.x() * multiplier);

    return mipmap[0](y,x) * multiplier * multiplier;
}


NORI_NAMESPACE_END
