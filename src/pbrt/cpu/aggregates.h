// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#ifndef PBRT_CPU_AGGREGATES_H
#define PBRT_CPU_AGGREGATES_H

#include <pbrt/pbrt.h>

#include <pbrt/cpu/primitive.h>
#include <pbrt/util/parallel.h>

#include <atomic>
#include <memory>
#include <vector>

namespace pbrt {

Primitive CreateAccelerator(const std::string &name, std::vector<Primitive> prims,
                            const ParameterDictionary &parameters);

struct BVHBuildNode;
struct BVHPrimitive;
struct LinearBVHNode;
struct MortonPrimitive;

// BVHAggregate Definition
class BVHAggregate {
  public:
    // BVHAggregate Public Types
    enum class SplitMethod { SAH, HLBVH, Middle, EqualCounts };

    // BVHAggregate Public Methods
    BVHAggregate(std::vector<Primitive> p, int maxPrimsInNode = 1,
                 SplitMethod splitMethod = SplitMethod::SAH);

    static BVHAggregate *Create(std::vector<Primitive> prims,
                                const ParameterDictionary &parameters);

    Bounds3f Bounds() const;
    pstd::optional<ShapeIntersection> Intersect(const Ray &ray, Float tMax) const;
    bool IntersectP(const Ray &ray, Float tMax) const;

  private:
    // BVHAggregate Private Methods
    BVHBuildNode *buildRecursive(ThreadLocal<Allocator> &threadAllocators,
                                 pstd::span<BVHPrimitive> bvhPrimitives,
                                 std::atomic<int> *totalNodes,
                                 std::atomic<int> *orderedPrimsOffset,
                                 std::vector<Primitive> &orderedPrims);
    BVHBuildNode *buildHLBVH(Allocator alloc,
                             const std::vector<BVHPrimitive> &primitiveInfo,
                             std::atomic<int> *totalNodes,
                             std::vector<Primitive> &orderedPrims);
    BVHBuildNode *emitLBVH(BVHBuildNode *&buildNodes,
                           const std::vector<BVHPrimitive> &primitiveInfo,
                           MortonPrimitive *mortonPrims, int nPrimitives, int *totalNodes,
                           std::vector<Primitive> &orderedPrims,
                           std::atomic<int> *orderedPrimsOffset, int bitIndex);
    BVHBuildNode *buildUpperSAH(Allocator alloc,
                                std::vector<BVHBuildNode *> &treeletRoots, int start,
                                int end, std::atomic<int> *totalNodes) const;
    int flattenBVH(BVHBuildNode *node, int *offset);

    // BVHAggregate Private Members
    int maxPrimsInNode;
    std::vector<Primitive> primitives;
    SplitMethod splitMethod;
    LinearBVHNode *nodes = nullptr;
};

struct KdTreeNode;
struct BoundEdge;

// KdTreeAggregate Definition
class KdTreeAggregate {
  public:
    // KdTreeAggregate Public Methods
    KdTreeAggregate(std::vector<Primitive> p, int isectCost = 5, int traversalCost = 1,
                    Float emptyBonus = 0.5, int maxPrims = 1, int maxDepth = -1);
    static KdTreeAggregate *Create(std::vector<Primitive> prims,
                                   const ParameterDictionary &parameters);
    pstd::optional<ShapeIntersection> Intersect(const Ray &ray, Float tMax) const;

    Bounds3f Bounds() const { return bounds; }

    bool IntersectP(const Ray &ray, Float tMax) const;

  private:
    // KdTreeAggregate Private Methods
    void buildTree(int nodeNum, const Bounds3f &bounds,
                   const std::vector<Bounds3f> &primBounds,
                   pstd::span<const int> primNums, int depth,
                   std::vector<BoundEdge> edges[3], pstd::span<int> prims0,
                   pstd::span<int> prims1, int badRefines);

    // KdTreeAggregate Private Members
    int isectCost, traversalCost, maxPrims;
    Float emptyBonus;
    std::vector<Primitive> primitives;
    std::vector<int> primitiveIndices;
    KdTreeNode *nodes;
    int nAllocedNodes, nextFreeNode;
    Bounds3f bounds;
};

struct GridVoxel {
    GridVoxel(std::vector<Primitive> &prims) : primitives(prims) {}

    void AddPrimitive(size_t prim_id) { prims_id.push_back(prim_id); }

    std::vector<size_t> prims_id;
    std::vector<Primitive> &primitives;
};

class GridAggregate {
  public:
    using VoxelCoordinate = Vector3i;
    // GridAggregate Public Methods
    GridAggregate(std::vector<Primitive> prims) : primitives(std::move(prims)) {
        calculateBounds();
        getVoxelSetting(bounds);
        createVoxels();
        size_t prim_id = 0;
        for (const auto &prim : primitives)
            fillVoxelExtents(prim, prim_id++);
    }

    static GridAggregate *Create(std::vector<Primitive> prims,
                                 const ParameterDictionary &parameters) {
        LOG_VERBOSE("GRID Creating!");
        LOG_VERBOSE("GRID Creating!");
        LOG_VERBOSE("GRID Creating!");
        return new GridAggregate(std::move(prims));
    }

    pstd::optional<ShapeIntersection> Intersect(const Ray &ray, Float tMax) const;

    Bounds3f Bounds() const { return bounds; }

    bool IntersectP(const Ray &ray, Float tMax) const {
        // TODO: Implement This
        for (const auto &p : primitives)
            if (p.IntersectP(ray, tMax))
                return true;
        return false;
    }

  private:
    // GridAggregate Private Methods
    void calculateBounds() {
        for (auto p : primitives)
            bounds = Union(bounds, p.Bounds());
    }
    // calculate required information for voxels
    void getVoxelSetting(const Bounds3f &bound) {
        // TODO
        delta = bounds.pMax - bounds.pMin;
        for (int axis = 0; axis < 3; ++axis) {
            // Hardcode to 64 voxels per axis
            nVoxel[axis] = 64;
            width[axis] = delta[axis] / nVoxel[axis];
            invWidth[axis] = 1.f / width[axis];
        }
    }
    // create Voxels for 3 dimension
    void createVoxels() {
        voxels.resize(nVoxel[0]);
        for (auto &x_voxels : voxels) {
            x_voxels.resize(nVoxel[1]);
            for (auto &xy_voxels : x_voxels)
                for (int i = 0; i < nVoxel[2]; ++i)
                    xy_voxels.emplace_back(primitives);
        }
    }

    void fillVoxelExtents(const Primitive &p, size_t prim_id) {
        auto primitive_bound = p.Bounds();
        auto voxel_min = posToVoxel(primitive_bound.pMin),
             voxel_max = posToVoxel(primitive_bound.pMax);
        for (auto x = voxel_min.x; x <= voxel_max.x; ++x)
            for (auto y = voxel_min.y; y <= voxel_max.y; ++y)
                for (auto z = voxel_min.z; z <= voxel_max.z; ++z) {
                    // // Debug
                    // std::string msg = "Add Primitive-" + std::to_string(prim_id) +
                    //                   " into voxel[" + std::to_string(x) + ", " +
                    //                   std::to_string(y) + ", " + std::to_string(z) +
                    //                   "]";
                    // LOG_VERBOSE(msg.c_str());
                    voxels[x][y][z].AddPrimitive(prim_id);
                }
    }

    VoxelCoordinate posToVoxel(const Point3f &p) {
        VoxelCoordinate voxel_coordinate;
        for (int axis = 0; axis < 3; ++axis)
            voxel_coordinate[axis] = posToVoxel(p, axis);
        return voxel_coordinate;
    }

    int posToVoxel(const Point3f &pos, int axis) const {
        auto v_pos = pstd::floor((pos[axis] - bounds.pMin[axis]) * invWidth[axis]);
        return Clamp(v_pos, 0, nVoxel[axis] - 1);
    }

    // GridAggregate Private Members
    std::vector<Primitive> primitives;

    std::vector<std::vector<std::vector<GridVoxel>>> voxels;

    Vector3f delta;
    Point3f nVoxel;
    Point3f width, invWidth;

    // bounds of all primitives
    Bounds3f bounds;
};

}  // namespace pbrt

#endif  // PBRT_CPU_AGGREGATES_H
