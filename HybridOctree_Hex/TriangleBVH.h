#ifndef TRIANGLE_BVH_H
#define TRIANGLE_BVH_H

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>
#include "Mesh.h"

// Broad phase only: callers retain their existing triangle predicates.
class TriangleBVH {
public:
    struct NearestResult {
        int triangle = -1;
        double distance = std::numeric_limits<double>::infinity();
        std::array<double, 3> point{{0, 0, 0}};
    };

    void Clear() { nodes.clear(); triangles.clear(); bounds.clear(); }
    bool Empty() const { return nodes.empty(); }

    void Build(const Mesh& mesh) {
        Clear();
        triangles.resize(mesh.eNum);
        std::iota(triangles.begin(), triangles.end(), 0);
        bounds.resize(mesh.eNum);
        for (int i = 0; i < mesh.eNum; ++i) {
            for (int axis = 0; axis < 3; ++axis) {
                double lo = mesh.v[mesh.e[i][0]][axis], hi = lo;
                for (int corner = 1; corner < 3; ++corner) {
                    lo = std::min(lo, mesh.v[mesh.e[i][corner]][axis]);
                    hi = std::max(hi, mesh.v[mesh.e[i][corner]][axis]);
                }
                // Conservative padding for the existing floating-point predicates.
                const double pad = 1e-10 * std::max(1.0, std::max(std::abs(lo), std::abs(hi)));
                bounds[i].lo[axis] = lo - pad;
                bounds[i].hi[axis] = hi + pad;
            }
        }
        nodes.reserve(mesh.eNum / 2 + 1);
        if (mesh.eNum) BuildNode(0, mesh.eNum);
    }

    // A finite segment, ray, or complete line, depending on [first, last].
    // Sorted IDs preserve the original narrow-phase traversal order.
    void LineCandidates(const double* point, const double* direction,
                        double first, double last, std::vector<int>& result) const {
        result.clear();
        if (!nodes.empty()) LineNode(0, point, direction, first, last, result);
        std::sort(result.begin(), result.end());
    }

    template<class Distance>
    NearestResult Nearest(const double* point, Distance distance) const {
        NearestResult result;
        if (!nodes.empty()) NearestNode(0, point, distance, result);
        return result;
    }

private:
    struct Box { std::array<double, 3> lo, hi; };
    struct Node {
        Box box;
        int first = 0, count = 0, left = -1, right = -1;
    };
    std::vector<Node> nodes;
    std::vector<int> triangles;
    std::vector<Box> bounds;

    int BuildNode(int first, int last) {
        const int index = static_cast<int>(nodes.size());
        nodes.emplace_back();
        Box box = bounds[triangles[first]];
        for (int i = first + 1; i < last; ++i)
            for (int axis = 0; axis < 3; ++axis) {
                box.lo[axis] = std::min(box.lo[axis], bounds[triangles[i]].lo[axis]);
                box.hi[axis] = std::max(box.hi[axis], bounds[triangles[i]].hi[axis]);
            }
        nodes[index].box = box;
        if (last - first <= 8) {
            nodes[index].first = first;
            nodes[index].count = last - first;
        } else {
            int axis = 0;
            for (int i = 1; i < 3; ++i)
                if (box.hi[i] - box.lo[i] > box.hi[axis] - box.lo[axis]) axis = i;
            const int middle = first + (last - first) / 2;
            std::nth_element(triangles.begin() + first, triangles.begin() + middle,
                             triangles.begin() + last, [&](int a, int b) {
                const double ca = bounds[a].lo[axis] + bounds[a].hi[axis];
                const double cb = bounds[b].lo[axis] + bounds[b].hi[axis];
                return ca < cb || (ca == cb && a < b);
            });
            const int left = BuildNode(first, middle);
            const int right = BuildNode(middle, last);
            nodes[index].left = left;
            nodes[index].right = right;
        }
        return index;
    }

    static bool HitsBox(const Box& box, const double* point, const double* direction,
                        double first, double last) {
        for (int axis = 0; axis < 3; ++axis) {
            if (direction[axis] == 0) {
                if (point[axis] < box.lo[axis] || point[axis] > box.hi[axis]) return false;
            } else {
                double a = (box.lo[axis] - point[axis]) / direction[axis];
                double b = (box.hi[axis] - point[axis]) / direction[axis];
                if (a > b) std::swap(a, b);
                first = std::max(first, a);
                last = std::min(last, b);
                if (first > last) return false;
            }
        }
        return true;
    }

    void LineNode(int index, const double* point, const double* direction,
                  double first, double last, std::vector<int>& result) const {
        const Node& node = nodes[index];
        if (!HitsBox(node.box, point, direction, first, last)) return;
        if (node.count) {
            for (int i = node.first; i < node.first + node.count; ++i) {
                const int triangle = triangles[i];
                if (HitsBox(bounds[triangle], point, direction, first, last))
                    result.push_back(triangle);
            }
        } else {
            LineNode(node.left, point, direction, first, last, result);
            LineNode(node.right, point, direction, first, last, result);
        }
    }

    static double BoxDistance(const Box& box, const double* point) {
        double squared = 0;
        for (int axis = 0; axis < 3; ++axis) {
            const double delta = std::max(0.0, std::max(box.lo[axis] - point[axis],
                                                       point[axis] - box.hi[axis]));
            squared += delta * delta;
        }
        return std::sqrt(squared);
    }

    template<class Distance>
    void NearestNode(int index, const double* point, Distance& distance,
                     NearestResult& result) const {
        const Node& node = nodes[index];
        if (BoxDistance(node.box, point) > result.distance) return;
        if (node.count) {
            for (int i = node.first; i < node.first + node.count; ++i) {
                const int triangle = triangles[i];
                if (BoxDistance(bounds[triangle], point) > result.distance) continue;
                std::array<double, 3> closest;
                const double d = distance(triangle, closest.data());
                // Lowest input ID wins exact ties, as in the original full scan.
                if (std::isfinite(d) && (d < result.distance ||
                    (d == result.distance && triangle < result.triangle))) {
                    result.distance = d;
                    result.triangle = triangle;
                    result.point = closest;
                }
            }
        } else {
            int first = node.left, second = node.right;
            if (BoxDistance(nodes[first].box, point) > BoxDistance(nodes[second].box, point))
                std::swap(first, second);
            NearestNode(first, point, distance, result);
            NearestNode(second, point, distance, result);
        }
    }
};
#endif
