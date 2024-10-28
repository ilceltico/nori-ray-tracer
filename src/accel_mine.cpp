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

#include <nori/accel.h>
#include <Eigen/Geometry>
#include <iostream>
#include <chrono>
#include <stdlib.h>
#include <list>
#include <unordered_map>
#include <tbb/task.h>
#include <tbb/task_scheduler_init.h>

NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh *mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

void Accel::build() {
    using std::chrono::high_resolution_clock;
    using std::chrono::duration;
    using std::chrono::milliseconds;


    // uint8_t maxDepth = 8;
    uint8_t maxDepth = (uint8_t) std::ceil((log(m_mesh->getTriangleCount()) / log(8))) + 1;
    
    #ifdef PARALLEL
    this->octreeSerialCutoff = (maxDepth + 1) / 2;
    #endif

    // system( "read -n 1 -s -p \"Press any key to continue...\"" ); //Debug
    auto t1 = high_resolution_clock::now();

    /* Create an initial list with all the triangle indices */
    std::vector<uint32_t> initialList = std::vector<uint32_t>();
    for (uint32_t i=0; i<m_mesh->getTriangleCount(); i++) {
        initialList.push_back(i);
    }

    #ifndef PARALLEL
    this->octreeRoot = buildOctreeSerial(m_bbox, initialList, maxDepth);
    #else
    tbb::task_scheduler_init init(OCTREE_THREAD_COUNT);
    this->octreeRoot = buildOctreeParallel(m_bbox, initialList, maxDepth);
    #endif

    std::vector<uint32_t>().swap(initialList); //Free memory

    #ifdef STATS
        maxDepthReached = maxDepth - maxDepthReached + 1;
    #endif

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> seconds = (t2 - t1) / 1000.0;

    /* Print stats */
    cout << "Octree build duration: " << seconds.count() << " seconds\n";
    cout << "Max depth set: " << (int) maxDepth << "\n";
    #ifdef PARALLEL
    cout << "Serial cutoff set: " << (int) this->octreeSerialCutoff << "\n";
    #endif
    #ifdef STATS
        cout << "Max depth reached: " << (int) maxDepthReached << "\n";
        cout << "Total number of nodes: " << leafNodes+parentNodes << "\n";
        cout << "Total number of parent nodes: " << parentNodes << "\n";
        cout << "Total number of leaf nodes: " << leafNodes << "\n";
        cout << "Total number of triangles in leaves " << trianglesInLeaves << "\n";
        cout << "Average number of triangles per leaf " << (double) trianglesInLeaves / leafNodes << "\n";
    #endif

    /* Debug */
    // OctreeNode *node = this->octreeRoot;
    // while(!node->isLeaf) {
    //     for (int i=0; i<8; i++) {
    //         if (((ParentOctreeNode*) node)->children[i]) {
    //             node = ((ParentOctreeNode*) node)->children[i];
    //             break;
    //         }
    //     }
    // }
    // LeafOctreeNode *leafNode = (LeafOctreeNode*) node;
    // cout << "Node size: " << sizeof(*leafNode) << "\n";
    // cout << "DEBUG: Number of triangles: " << ((LeafOctreeNode*) node)->triangleIndices.size() << "\n";
    // cout << "DEBUG: Number of triangles: " << ((LeafOctreeNode*) node)->boundingBox.toString() << "\n";


    // system( "read -n 1 -s -p \"Press any key to continue...\"" );

    // ParentOctreeNode *parentNode = (ParentOctreeNode*) this->octreeRoot;
    // cout << "Node size: " << sizeof(*parentNode) << "\n";

}

#ifdef STATS
Accel::OctreeNode * Accel::buildOctreeSerial(BoundingBox3f const &boundingBox, std::vector<uint32_t> &triangleIndices, uint8_t maxDepth) {
#else
Accel::OctreeNode * Accel::buildOctreeSerial(BoundingBox3f const &boundingBox, std::vector<uint32_t> &triangleIndices, uint8_t maxDepth) const {
#endif 
    // cout << "Depth: " << depth << "\n"; //Debug

    /* No indices, empty node */
    if (triangleIndices.empty()) {
        return nullptr;
    }

    /* Leaf node */
    if (triangleIndices.size() < OCTREE_MAX_TRIANGLES || maxDepth <= (uint8_t) 1) {
        #ifdef STATS
            leafNodes++;
            trianglesInLeaves += triangleIndices.size();
            if (this->maxDepthReached > maxDepth) {
                this->maxDepthReached = maxDepth;
            }
        #endif
        // cout << "Number of leaf nodes: " << leafNodes << ", Number of leaf triangles: " << trianglesInLeaves << "\n"; //Debug
        return new LeafOctreeNode(boundingBox, triangleIndices);
    }

    /* Create and populate triangle lists for all the children */
    std::vector<uint32_t>* childrenTriangleLists[8];
    BoundingBox3f* childrenBoundingBoxes[8];
    Point3f center = boundingBox.getCenter();
    for (int i=0; i<8; i++) {
        childrenTriangleLists[i] = new std::vector<uint32_t>();
        childrenBoundingBoxes[i] = new BoundingBox3f(center);
        childrenBoundingBoxes[i]->expandBy(boundingBox.getCorner(i));
    }
    for (auto it=triangleIndices.begin(); it != triangleIndices.end(); ++it) {
        //Get the triangle Bounding Box
        BoundingBox3f triangleBoundingBox = m_mesh->getBoundingBox(*it);
        for (int i=0; i<8; i++) {
            if (triangleBoundingBox.overlaps(*childrenBoundingBoxes[i])) {
                childrenTriangleLists[i]->push_back(*it);
            }
        }
    }

    /* Recursion on child nodes */
    ParentOctreeNode *node = new ParentOctreeNode(boundingBox);
    for (int i=0; i<8; i++) {
        node->children[i] = buildOctreeSerial(*childrenBoundingBoxes[i], *childrenTriangleLists[i], maxDepth-1);

        //They are copied in the constructor of an OctreeNode. Avoiding the copy does improve times a little, but makes memory management more difficult.
        delete childrenTriangleLists[i];
        delete childrenBoundingBoxes[i];
    }

    #ifdef STATS
        parentNodes++;
    #endif
    return node;
}


#ifdef PARALLEL
Accel::OctreeNode* Accel::buildOctreeParallel(BoundingBox3f const &boundingBox, std::vector<uint32_t> &triangleIndices, uint8_t maxDepth) const {
    OctreeNode* result;

    BuildOctreeTask& task = *new(tbb::task::allocate_root()) BuildOctreeTask(this, boundingBox, triangleIndices, maxDepth, result);
    tbb::task::spawn_root_and_wait(task);

    return result;
}

tbb::task* Accel::BuildOctreeTask::execute() {

    /* Cutoff threashold reached, execute till maxDepth serially */
    if (maxDepth <= accel->octreeSerialCutoff) {
        resultNode = accel->buildOctreeSerial(boundingBox, triangleIndices, maxDepth);
    }

    else {

         /* No indices, empty node */
        if (triangleIndices.empty()) {
            return NULL;
        }

        /* Leaf node */
        if (triangleIndices.size() <= OCTREE_MAX_TRIANGLES || maxDepth <= (uint8_t) 1) {
            resultNode = new LeafOctreeNode(boundingBox, triangleIndices);
        }

        /* Create and populate triangle lists for all the children */
        std::vector<uint32_t>* childrenTriangleLists[8];
        BoundingBox3f* childrenBoundingBoxes[8];
        Point3f center = boundingBox.getCenter();
        for (int i=0; i<8; i++) {
            childrenTriangleLists[i] = new std::vector<uint32_t>();
            childrenBoundingBoxes[i] = new BoundingBox3f(center);
            childrenBoundingBoxes[i]->expandBy(boundingBox.getCorner(i));
        }
        for (auto it=triangleIndices.begin(); it != triangleIndices.end(); ++it) {
            //Get the triangle Bounding Box
            BoundingBox3f triangleBoundingBox = accel->m_mesh->getBoundingBox(*it);
            for (int i=0; i<8; i++) {
                if (triangleBoundingBox.overlaps(*childrenBoundingBoxes[i])) {
                    childrenTriangleLists[i]->push_back(*it);
                }
            }
        }

        /* Recursion on child nodes */
        resultNode = new ParentOctreeNode(boundingBox);
        BuildOctreeTask* tasks[8];
        for (uint8_t i=0; i<8; i++) {
            tasks[i] = new (allocate_child()) BuildOctreeTask(this->accel, *childrenBoundingBoxes[i], *childrenTriangleLists[i], maxDepth-1, ((ParentOctreeNode*) resultNode)->children[i]);
            
            //They are copied in the constructor of a BuildOctreeTask. Avoiding the copy does improve times a little, but makes memory management more difficult.
            delete childrenTriangleLists[i];
            delete childrenBoundingBoxes[i];
        }
        set_ref_count(9);
        for (uint8_t i=0; i<7; i++) {
            spawn(*tasks[i]);
        }
        /* Wait for all the children tasks to return */
        spawn_and_wait_for_all(*tasks[7]);
    }
    return NULL;
}
#endif


bool Accel::findIntersection(Ray3f &ray, Intersection &its, bool shadowRay, OctreeNode *node, uint32_t &f, float &u, float &v, float &t, float &nearT, float &farT) const {
    bool foundIntersection = false;

    /* Leaf node */
    if (node->isLeaf) {            
        for (auto const& it : ((LeafOctreeNode*) node)->triangleIndices) {
            if (m_mesh->rayIntersect(it, ray, u, v, t)) {
                /* An intersection was found! Can terminate
                immediately if this is a shadow ray query */
                if (shadowRay)
                    return true;
                ray.maxt = its.t = t;
                its.uv = Point2f(u, v);
                its.mesh = m_mesh;
                f = it;
                foundIntersection = true;
            }
        }

    /* Parent node */
    } else {                      
        uint8_t i, j = 0;
        uint8_t childrenIndicesToQueue[8] = {0,1,2,3,4,5,6,7};
        float nearTmap[8] = {0,0,0,0,0,0,0,0};

        /* Check which children can be intersected by the current ray and cache the nearT information */
        for (i=0; i<8; i++) {
            if (((ParentOctreeNode*) node)->children[i] && ((ParentOctreeNode*) node)->children[i]->boundingBox.rayIntersect(ray, nearT, farT)) {
                if (ray.mint <= farT && nearT <= ray.maxt) {
                    childrenIndicesToQueue[j++] = i;
                    nearTmap[i] = nearT;
                }
            }
        }

        /* Sort the children by nearT using the cached values */
        std::sort(childrenIndicesToQueue, childrenIndicesToQueue+j, [&nearTmap](uint8_t a, uint8_t b) {
                                                                                    return nearTmap[a] < nearTmap[b];
                                                                                });

        /* Recursion on the sorted list of children */
        for (i=0; i<j; i++) {
            if (findIntersection(ray, its, shadowRay, ((ParentOctreeNode*) node)->children[childrenIndicesToQueue[i]], f, u, v, t, nearT, farT)) {
                foundIntersection = true;
                if (shadowRay)
                    return true;
            }
            /* Avoiding recursion if, after exploring the previous children, the following children can be pruned */
            if (i<j-1 && nearTmap[childrenIndicesToQueue[i+1]] > ray.maxt)
                break;
        }
    }
    return foundIntersection;
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    /* Brute force search through all triangles */
    // for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
    //     float u, v, t;
    //     if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
    //         /* An intersection was found! Can terminate
    //            immediately if this is a shadow ray query */
    //         if (shadowRay)
    //             return true;
    //         ray.maxt = its.t = t;
    //         its.uv = Point2f(u, v);
    //         its.mesh = m_mesh;
    //         f = idx;
    //         foundIntersection = true;
    //     }
    // }

    /* Iterative Depth-first search*/
    // std::list<OctreeNode*> queue = {this->octreeRoot};
    // // std::list<BoundingBox3f> queueBB = {this->m_bbox};
    // OctreeNode* currentNode;
    // // BoundingBox3f currentBB;
    // uint32_t idx;
    // // int leaves = 0;
    // while (!queue.empty()) {
    //     currentNode = queue.front();
    //     queue.pop_front();
    //     // currentBB = queueBB.front();
    //     // queueBB.pop_front();

    //     // if (!currentNode) {
    //     //     continue;
    //     // }

    //     // if (!currentBB.rayIntersect(ray)) {
    //     //     continue;
    //     // }

    //     if (currentNode->isLeaf()) {    //Leaf node
    //         // leaves++;
    //         for (uint32_t i = 0; i < ((LeafOctreeNode*) currentNode)->triangleIndices.size(); ++i) {
    //             idx = ((LeafOctreeNode*) currentNode)->triangleIndices[i];
    //             float u, v, t;
    //             if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
    //                 /* An intersection was found! Can terminate
    //                 immediately if this is a shadow ray query */
    //                 if (shadowRay)
    //                     return true;
    //                 ray.maxt = its.t = t;
    //                 its.uv = Point2f(u, v);
    //                 its.mesh = m_mesh;
    //                 f = idx;
    //                 foundIntersection = true;
    //             }
    //         }
    //     } else {                        //Parent node
    //         for (int i=0; i<8; i++) {
    //             // BoundingBox3f nextBB = BoundingBox3f(currentBB.getCenter());
    //             // nextBB.expandBy(currentBB.getCorner(i));
    //             // if (((ParentOctreeNode*) currentNode)->children[i] && nextBB.rayIntersect(ray)) {
    //             if (((ParentOctreeNode*) currentNode)->children[i] && ((ParentOctreeNode*) currentNode)->children[i]->boundingBox.rayIntersect(ray)) {
    //                 queue.push_front(((ParentOctreeNode*) currentNode)->children[i]);
    //                 // queueBB.push_front(nextBB);
    //             }
    //         }
    //     }
    // }


    /* Iterative Depth-first search with distance heuristic*/
    // float nearT;
    // float farT;
    // OctreeNodeWithDistance currentNode = OctreeNodeWithDistance(this->octreeRoot, this->m_bbox.rayIntersect(ray, nearT, farT));
    // std::list<OctreeNodeWithDistance> queue = {currentNode};
    // ParentOctreeNode* parentCurrentNode;
    // LeafOctreeNode* leafCurrentNode;
    // uint32_t idx, i, j;
    // float u, v, t;
    // std::vector<uint8_t> childrenIndicesToQueue = {0,1,2,3,4,5,6,7};
    // float nearTmap[8];
    // // int leaves = 0;
    // while (!queue.empty()) {
    //     currentNode = queue.front();
    //     queue.pop_front();

    //     if (currentNode.nearT > ray.maxt) {
    //         continue;
    //     }

    //     // currentBB = queueBB.front();
    //     // queueBB.pop_front();

    //     // if (!currentNode) {
    //     //     continue;
    //     // }

    //     // if (!currentBB.rayIntersect(ray)) {
    //     //     continue;
    //     // }

    //     if (currentNode.node->isLeaf()) {    //Leaf node
    //         // leaves++;
    //         leafCurrentNode = ((LeafOctreeNode*) currentNode.node);
    //         for (i = 0; i < leafCurrentNode->triangleIndices.size(); ++i) {
    //             idx = leafCurrentNode->triangleIndices[i];
    //             if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
    //                 /* An intersection was found! Can terminate
    //                 immediately if this is a shadow ray query */
    //                 if (shadowRay)
    //                     return true;
    //                 ray.maxt = its.t = t;
    //                 its.uv = Point2f(u, v);
    //                 its.mesh = m_mesh;
    //                 f = idx;
    //                 foundIntersection = true;
    //             }
    //         }
    //     } else {                        //Parent node
    //         // std::vector<int> childrenIndicesToQueue;
    //         // std::unordered_map<int, float> nearTmap;
    //         // for (int i=0; i<8; i++) {
    //         //     if (((ParentOctreeNode*) currentNode)->children[i] && ((ParentOctreeNode*) currentNode)->children[i]->boundingBox.rayIntersect(ray, nearT, farT)) {
    //         //         if (ray.mint <= farT && nearT <= ray.maxt) {
    //         //             childrenIndicesToQueue.push_back(i);
    //         //             nearTmap.insert(std::pair<int, float>(i,nearT));
    //         //         }
    //         //     }
    //         // }
    //         // std::sort(childrenIndicesToQueue.begin(), childrenIndicesToQueue.end(), [&nearTmap](int a, int b) {
    //         //                                                                             return nearTmap.find(a)->second > nearTmap.find(b)->second;
    //         //                                                                         });
    //         // for (int childIndex : childrenIndicesToQueue) {
    //         //     queue.push_front(((ParentOctreeNode*) currentNode)->children[childIndex]);
    //         // }

    //         parentCurrentNode = ((ParentOctreeNode*) currentNode.node);
    //         j = 0;
    //         for (i=0; i<8; i++) {
    //             if (parentCurrentNode->children[i] && parentCurrentNode->children[i]->boundingBox.rayIntersect(ray, nearT, farT)) {
    //                 if (ray.mint <= farT && nearT <= ray.maxt) {
    //                     childrenIndicesToQueue[j++] = i;
    //                     nearTmap[i] = nearT;
    //                 }
    //             }
    //         }
    //         std::sort(childrenIndicesToQueue.begin(), childrenIndicesToQueue.begin()+j, [&nearTmap](uint8_t a, uint8_t b) {
    //                                                                                     return nearTmap[a] > nearTmap[b];
    //                                                                                 });
    //         for (i=0; i<j; i++) {
    //             queue.push_front(OctreeNodeWithDistance(parentCurrentNode->children[childrenIndicesToQueue[i]], nearTmap[childrenIndicesToQueue[i]]));
    //         }
    //     }
    // }

    /* Recursive Depth-first */
    float u, v, t;
    float nearT, farT;
    foundIntersection = findIntersection(ray, its, shadowRay, this->octreeRoot, f, u, v, t, nearT, farT);

    if (shadowRay)
        return foundIntersection;

    // cout << "Leaves hit by the current ray: " << leaves << "\n";
    // system( "read -n 1 -s -p \"Press any key to continue...\"" );


    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                 bary.y() * N.col(idx1) +
                 bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

NORI_NAMESPACE_END

