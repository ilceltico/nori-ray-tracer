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

#include <nori/mesh.h>
#include <nori/bbox.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/warp.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

Mesh::Mesh() { }

Mesh::~Mesh() {
    delete m_bsdf;
    delete m_emitter;
}

void Mesh::activate() {
    if (!m_bsdf) {
        /* If no material was assigned, instantiate a diffuse BRDF */
        m_bsdf = static_cast<BSDF *>(
            NoriObjectFactory::createInstance("diffuse", PropertyList()));
        m_bsdf->activate();
    }

    /* Build a PDF to sample triangles based on their surface area */
    m_trianglePDF = DiscretePDF(getTriangleCount());
    for (uint32_t i=0; i<getTriangleCount(); i++) {
        m_trianglePDF.append(surfaceArea(i));
    }
    totalSurfaceArea = m_trianglePDF.normalize();
    inverseTotalSurfaceArea = 1.0f / totalSurfaceArea;

    /* Initialize the Emitter, if any */
    if (m_emitter) {
        if (m_emitter->isAssociatedToMesh()) {
            m_emitter->setAssociatedMesh(this);
        }
    }

    /* Compute and store the dpdu and dpdv vectors, needed for some BSDFs */
    
    if (m_bsdf->requiresDifferentials()) {

        m_dpdu.resize(3, m_F.cols());
        m_dpdv.resize(3, m_F.cols());

        for (long i=0; i<m_F.cols(); i++) {
            uint32_t idx0 = m_F(0,i),
                    idx1 = m_F(1,i),
                    idx2 = m_F(2,i);

            const Point3f
                &v0 = m_V.col(idx0),
                &v1 = m_V.col(idx1),
                &v2 = m_V.col(idx2);

            const Point2f
                &uv0 = m_UV.col(idx0),
                &uv1 = m_UV.col(idx1),
                &uv2 = m_UV.col(idx2);

            Vector3f dP1 = v1 - v0, dP2 = v2 - v0;
            Vector2f dUV1 = uv1 - uv0, dUV2 = uv2 - uv0;
            Normal3f n = Normal3f(dP1.cross(dP2));
            float length = sqrt(n.x()*n.x() + n.y()*n.y() + n.z()*n.z());
            if (length == 0) {
                continue;
            }

            float determinant = dUV1.x() * dUV2.y() - dUV1.y() * dUV2.x();
            if (determinant == 0) {
                /* The user-specified parameterization is degenerate. Pick
                arbitrary tangents that are perpendicular to the geometric normal */
                Vector3f dpdu, dpdv;
                coordinateSystem(n.normalized(), dpdu, dpdv);
                m_dpdu.col(i) = dpdu;
                m_dpdv.col(i) = dpdv;
            } else {
                float invDet = 1.0f / determinant;
                m_dpdu.col(i) = ( dUV2.y() * dP1 - dUV1.y() * dP2) * invDet;
                m_dpdv.col(i) = (-dUV2.x() * dP1 + dUV1.x() * dP2) * invDet;
            }
        }
    }
   
}

float Mesh::surfaceArea(uint32_t index) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);

    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    return 0.5f * Vector3f((p1 - p0).cross(p2 - p0)).norm();
}

bool Mesh::rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t, Vector3f &dpdu, Vector3f &dpdv) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    /* Find vectors for two edges sharing v[0] */
    Vector3f edge1 = p1 - p0, edge2 = p2 - p0;

    /* Begin calculating determinant - also used to calculate U parameter */
    Vector3f pvec = ray.d.cross(edge2);

    /* If determinant is near zero, ray lies in plane of triangle */
    float det = edge1.dot(pvec);

    if (det > -1e-8f && det < 1e-8f)
        return false;
    float inv_det = 1.0f / det;

    /* Calculate distance from v[0] to ray origin */
    Vector3f tvec = ray.o - p0;

    /* Calculate U parameter and test bounds */
    u = tvec.dot(pvec) * inv_det;
    if (u < 0.0 || u > 1.0)
        return false;

    /* Prepare to test V parameter */
    Vector3f qvec = tvec.cross(edge1);

    /* Calculate V parameter and test bounds */
    v = ray.d.dot(qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0)
        return false;

    /* Ray intersects triangle -> compute t */
    t = edge2.dot(qvec) * inv_det;

    if (m_bsdf->requiresDifferentials()) {
        dpdu = m_dpdu.col(index);
        dpdv = m_dpdv.col(index);
    }

    return t >= ray.mint && t <= ray.maxt;
}

BoundingBox3f Mesh::getBoundingBox(uint32_t index) const {
    BoundingBox3f result(m_V.col(m_F(0, index)));
    result.expandBy(m_V.col(m_F(1, index)));
    result.expandBy(m_V.col(m_F(2, index)));
    return result;
}

Point3f Mesh::getCentroid(uint32_t index) const {
    return (1.0f / 3.0f) *
        (m_V.col(m_F(0, index)) +
         m_V.col(m_F(1, index)) +
         m_V.col(m_F(2, index)));
}

void Mesh::addChild(NoriObject *obj) {
    switch (obj->getClassType()) {
        case EBSDF:
            if (m_bsdf)
                throw NoriException(
                    "Mesh: tried to register multiple BSDF instances!");
            m_bsdf = static_cast<BSDF *>(obj);
            break;

        case EEmitter: {
                Emitter *emitter = static_cast<Emitter *>(obj);
                if (m_emitter)
                    throw NoriException(
                        "Mesh: tried to register multiple Emitter instances!");
                m_emitter = emitter;
            }
            break;

        case EMedium: {
                Medium *medium = static_cast<Medium *>(obj);
                if (medium->isInternalMedium()) {
                    if (m_medium_internal)
                        throw NoriException(
                            "Mesh: tried to register multiple internal Medium instances!");
                    m_medium_internal = static_cast<Medium *>(obj);
                } else {
                    if (m_medium_external)
                        throw NoriException(
                            "Mesh: tried to register multiple external Medium instances!");
                    m_medium_external = static_cast<Medium *>(obj);
                }
            }
            break;

        default:
            throw NoriException("Mesh::addChild(<%s>) is not supported!",
                                classTypeName(obj->getClassType()));
    }
}

std::string Mesh::toString() const {
    return tfm::format(
        "Mesh[\n"
        "  name = \"%s\",\n"
        "  vertexCount = %i,\n"
        "  triangleCount = %i,\n"
        "  bsdf = %s,\n"
        "  emitter = %s\n"
        "  internal medium = %s\n"
        "  external medium = %s\n"
        "]",
        m_name,
        m_V.cols(),
        m_F.cols(),
        m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
        m_emitter ? indent(m_emitter->toString()) : std::string("null"),
        m_medium_internal ? indent(m_medium_internal->toString()) : std::string("null"),
        m_medium_external ? indent(m_medium_external->toString()) : std::string("null")
    );
}

float Mesh::uniformSquareToUniformSurface(const Point2f &sample, Point3f &surfacePoint, Normal3f &surfaceNormal, float &u, float &v) const {
    /* Sample a triangle based on surface area */
    float x = sample.x();
    uint32_t triangleIndex = m_trianglePDF.sampleReuse(x);

    /* Get barycentric coefficients from uniform 2D sample */
    float sqrt1x = sqrt(1 - x);
    float alpha = 1 - sqrt1x;
    float beta = sample.y() * sqrt1x;

    uint32_t vertexIndex1, vertexIndex2, vertexIndex3;
    vertexIndex1 = m_F(0, triangleIndex);
    vertexIndex2 = m_F(1, triangleIndex);
    vertexIndex3 = m_F(2, triangleIndex);

    const Point3f p1 = m_V.col(vertexIndex1), p2 = m_V.col(vertexIndex2), p3 = m_V.col(vertexIndex3);

    /* Get 3D coordinates from the barycentric coefficients */
    surfacePoint = alpha * p1 + beta * p2 + (1 - alpha - beta) * p3;

    Vector3f edge1 = p2 - p1;
    Vector3f edge2 = p3 - p1;

    /* Get surface normal on sampled point */
    /* If the vertex normals are available, interpolate them */
    if (m_N.size() > 0) {
        surfaceNormal = alpha * m_N.col(vertexIndex1) + beta * m_N.col(vertexIndex2) + (1 - alpha - beta) * m_N.col(vertexIndex3);
    }
    /* Else, use the face normal */
    else {
        surfaceNormal = edge1.cross(edge2);
    }
    surfaceNormal.normalize();




    // Compute also the UV coordinates of the sampled point, in case of textured lights

    /* Begin calculating determinant - also used to calculate U parameter */
    Vector3f pvec = (-surfaceNormal).cross(edge2);

    /* If determinant is near zero, ray lies in plane of triangle */
    float det = edge1.dot(pvec);

    if (det > -1e-8f && det < 1e-8f)
        return false;
    float inv_det = 1.0f / det;

    /* Calculate distance from v[0] to ray origin */
    Vector3f tvec = (surfacePoint - surfaceNormal) - p1;

    /* Calculate U parameter and test bounds */
    u = tvec.dot(pvec) * inv_det;

    /* Prepare to test V parameter */
    Vector3f qvec = tvec.cross(edge1);

    /* Calculate V parameter and test bounds */
    v = (-surfaceNormal).dot(qvec) * inv_det;


    return inverseTotalSurfaceArea;
}


std::string Intersection::toString() const {
    if (!mesh)
        return "Intersection[invalid]";

    return tfm::format(
        "Intersection[\n"
        "  p = %s,\n"
        "  t = %f,\n"
        "  uv = %s,\n"
        "  shFrame = %s,\n"
        "  geoFrame = %s,\n"
        "  mesh = %s\n"
        "]",
        p.toString(),
        t,
        uv.toString(),
        indent(shFrame.toString()),
        indent(geoFrame.toString()),
        mesh ? mesh->toString() : std::string("null")
    );
}

NORI_NAMESPACE_END
