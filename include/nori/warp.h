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

#pragma once

#include <nori/common.h>
#include <nori/sampler.h>
#include <nori/bitmap.h>

NORI_NAMESPACE_BEGIN

class Mipmap {
public:
    Mipmap() {}

    Mipmap(std::string filename) {
        buildMipmap(filename);
    }

    void buildMipmap(std::string filename);

    Point2f squareToHierarchicalSampleWarping(const Point2f &sample) const;

    float squareToHierarchicalSampleWarpingPdf(const Point2f &p) const;

    int getCols() const {
        return mipmap[0].cols();
    }

    int getRows() const {
        return mipmap[0].rows();
    }

private:
    std::vector<MatrixXf> mipmap;
};


/// A collection of useful warping functions for importance sampling
class Warp {
public:
    /// Dummy warping function: takes uniformly distributed points in a square and just returns them
    static Point2f squareToUniformSquare(const Point2f &sample);

    /// Probability density of \ref squareToUniformSquare()
    static float squareToUniformSquarePdf(const Point2f &p);

    /// Sample a 2D tent distribution
    static Point2f squareToTent(const Point2f &sample);

    /// Probability density of \ref squareToTent()
    static float squareToTentPdf(const Point2f &p);

    /// Uniformly sample a vector on a 2D disk with radius 1, centered around the origin
    static Point2f squareToUniformDisk(const Point2f &sample);

    /// Probability density of \ref squareToUniformDisk()
    static float squareToUniformDiskPdf(const Point2f &p);

    /// Uniformly sample a vector on the unit sphere with respect to solid angles
    static Vector3f squareToUniformSphere(const Point2f &sample);

    /// Probability density of \ref squareToUniformSphere()
    static float squareToUniformSpherePdf(const Vector3f &v);

    /// Uniformly sample a vector on the unit hemisphere around the pole (0,0,1) with respect to solid angles
    static Vector3f squareToUniformHemisphere(const Point2f &sample);

    /// Probability density of \ref squareToUniformHemisphere()
    static float squareToUniformHemispherePdf(const Vector3f &v);

    /// Uniformly sample a vector on the unit hemisphere around the pole (0,0,1) with respect to projected solid angles
    static Vector3f squareToCosineHemisphere(const Point2f &sample);

    /// Probability density of \ref squareToCosineHemisphere()
    static float squareToCosineHemispherePdf(const Vector3f &v);

    /// Warp a uniformly distributed square sample to a Beckmann distribution * cosine for the given 'alpha' parameter
    static Vector3f squareToBeckmann(const Point2f &sample, float alpha);

    /// Probability density of \ref squareToBeckmann()
    static float squareToBeckmannPdf(const Vector3f &m, float alpha);

    /// Warp a uniformly distributed square sample using Hierarchical Sample Warping to an image defining the sampling density
    static Point2f squareToHierarchicalSampleWarping(const Point2f &sample, const Mipmap &mipmap);

    /// Probability density of \ref squareToHierarchicalSampleWarping()
    static float squareToHierarchicalSampleWarpingPdf(const Point2f &p, const Mipmap &mipmap);

};

NORI_NAMESPACE_END
