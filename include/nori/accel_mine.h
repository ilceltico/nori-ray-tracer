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

#include <nori/mesh.h>
#include <vector>
#include <list>
#include <tbb/task.h>

// #define STATS //Use only *without* the PARALLEL flag. Updates and prints Octree stats in a non concurrency-friendly way.
#define PARALLEL //Sets Octree building to multi-threaded using TBB

NORI_NAMESPACE_BEGIN

const uint32_t OCTREE_MAX_TRIANGLES = 10;
#ifdef PARALLEL
const uint8_t OCTREE_THREAD_COUNT = 8;
#endif

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */
class Accel {
public:
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh *mesh);

    /// Build the acceleration data structure (currently a no-op)
    void build();

    /// Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const { return m_bbox; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its, bool shadowRay) const;

private:

    /**
     * \brief Abstract struct to be used as a base for Octree nodes.
     */
    struct OctreeNode{
        BoundingBox3f boundingBox; //24 Bytes
        bool isLeaf = true; //Most nodes are leaves
        // virtual bool isLeaf() = 0;
    };

    /**
     * \brief A Parent Octree node has exactly 8 children, and \c isLeaf=false .
     */
    struct ParentOctreeNode : public OctreeNode {
        OctreeNode* children[8];    //64 Bytes

        ParentOctreeNode(BoundingBox3f boundingBox) : children() {
            this->boundingBox = boundingBox;
            this->isLeaf = false;
        }

        // bool isLeaf() {             //8 Bytes
        //     return false;
        // }
    };

    /**
     * \brief A Leaf node does not have any child, but it has a vector of triangles
     * (expressed as integer indices for the Mesh). \c isLeaf=true .
     */
    struct LeafOctreeNode : public OctreeNode {
        std::vector<uint32_t> triangleIndices; //24 Bytes

        LeafOctreeNode(BoundingBox3f boundingBox, std::vector<uint32_t> &triangleIndices) {
            this->triangleIndices = triangleIndices; //Copy
            this->boundingBox = boundingBox;
        }

        // bool isLeaf() {                           //8 Bytes
        //     return true;
        // }
    };


    #ifdef PARALLEL
    /**
     * \brief Helper class for the parallel construction of the Octree
    */
    class BuildOctreeTask : public tbb::task {
        public:
            Accel const *accel;
            BoundingBox3f boundingBox;
            std::vector<uint32_t> triangleIndices;
            uint8_t maxDepth;

            OctreeNode* &resultNode;

            BuildOctreeTask(Accel const *accel, BoundingBox3f boundingBox, std::vector<uint32_t> &triangleIndices, uint8_t maxDepth, OctreeNode* &resultNode_) : resultNode(resultNode_) {
                this->accel = accel;
                this->boundingBox = boundingBox;
                this->triangleIndices = triangleIndices;
                this->maxDepth = maxDepth;
            }

            task* execute();
    };
    #endif


    
    #ifdef STATS
    /**
     * \brief Private method to be used recursively for the single-threaded construction of the Octree
     * 
     * \param boundingBox
     *      The BoundingBox for the currently evaluated node. To be set as the whole mesh's BoundingBox on first call.
     * 
     * \param triangleList
     *      Vector of triangle indices in the current node. To be set as the whole vector of the mesh's triangle indices on first call.
     * 
     * \param maxDepth
     *      Maximum depth allowed for the Octree. It is decreased at each recursion step. A \c maxDepth=1 means that only the root node is created.
     * 
     */
    OctreeNode* buildOctreeSerial(BoundingBox3f const &boundingBox, std::vector<uint32_t> &triangleList, uint8_t maxDepth);
    #else
    /**
     * \brief Private method to be used recursively for the single-threaded construction of the Octree
     * 
     * \param boundingBox
     *      The BoundingBox for the currently evaluated node. To be set as the whole mesh's BoundingBox on first call.
     * 
     * \param triangleList
     *      Vector of triangle indices in the current node. To be set as the whole vector of the mesh's triangle indices on first call.
     * 
     * \param maxDepth
     *      Maximum depth allowed for the Octree. It is decreased at each recursion step. A \c maxDepth=1 means that only the root node is created.
     * 
     */
    OctreeNode* buildOctreeSerial(BoundingBox3f const &boundingBox, std::vector<uint32_t> &triangleList, uint8_t maxDepth) const;
    #endif

    #ifdef PARALLEL
    /**
     * \brief Private method to be for a TBB-based multi-threaded construction of the Octree
     * 
     * \param boundingBox
     *      The BoundingBox for the currently evaluated node. To be set as the whole mesh's BoundingBox on first call.
     * 
     * \param triangleList
     *      Vector of triangle indices in the current node. To be set as the whole vector of the mesh's triangle indices on first call.
     * 
     * \param maxDepth
     *      Maximum depth allowed for the Octree. It is decreased at each recursion step. A \c maxDepth=1 means that only the root node is created.
     * 
     */
    OctreeNode* buildOctreeParallel(BoundingBox3f const &boundingBox, std::vector<uint32_t> &triangleList, uint8_t maxDepth) const;
    #endif

    /**
     * \brief Private method to be used recursively for the octree-based ray traversal
     * 
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     * 
     * \param node
     *    The OctreeNode to be evaluated in this call
     * 
     * \param f
     *    At finished execution, it will contain the triangle index of the closest intersection
     * 
     * \param u
     *    At finished execution, it will contain the 'U' component of the intersection in barycentric coordinates
     * 
     * \param v
     *    At finished execution, it will contain the 'V' component of the intersection in barycentric coordinates
     * 
     * \param t
     *    At finished execution, it will contain the distance from the ray origin to the intersection point
     * 
     * \param nearT
     *    Parameter used in the function to avoid unnecessary memory instantiation
     * 
     * \param farT
     *    Parameter used in the function to avoid unnecessary memory instantiation
     *
     * \return \c true if an intersection was found
     * 
    */
    bool findIntersection(Ray3f &ray, Intersection &its, bool shadowRay, OctreeNode *node, uint32_t &f, float &u, float &v, float &t, float &nearT, float &farT) const;


    Mesh         *m_mesh = nullptr; ///< Mesh (only a single one for now)
    BoundingBox3f m_bbox;           ///< Bounding box of the entire scene

    OctreeNode   *octreeRoot = nullptr; ///< Root node of the Octree for faster ray operations

    #ifdef PARALLEL
    uint8_t octreeSerialCutoff;     /// Depth at which executing is reverted to serial
    #endif

    /// Stats variables 
    #ifdef STATS
        uint32_t      leafNodes = 0;
        uint32_t      parentNodes = 0;
        uint32_t      trianglesInLeaves = 0;
        uint8_t       maxDepthReached = -1;
    #endif

};

NORI_NAMESPACE_END
