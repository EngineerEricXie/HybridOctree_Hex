#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
// Test the existing private stages without widening the application API.
#define private public
#include "../HexGen.cpp"
#undef private

static void Require(bool value, const std::string& message) {
    if (!value) throw std::runtime_error(message);
}

static void CubeSurface(hexGen& gen) {
    gen.triMesh.Initialize(12, 8, 3);
    gen.triMesh.vNum = 8;
    gen.triMesh.eNum = 12;
    const double vertices[8][3] = {{0,0,0}, {100,0,0}, {100,100,0}, {0,100,0},
                                  {0,0,100}, {100,0,100}, {100,100,100}, {0,100,100}};
    std::memcpy(gen.triMesh.v, vertices, sizeof(vertices));
    for (int face = 0; face < 6; ++face) {
        const int indices[2][3] = {{0,1,2}, {0,2,3}};
        for (int t = 0; t < 2; ++t)
            for (int j = 0; j < 3; ++j)
                gen.triMesh.e[2 * face + t][j] = fIdC[face][indices[t][j]];
    }
    gen.EnsureSurfaceIndex();
}

static void Geometry(hexGen& gen) {
    gen.EnsureSurfaceIndex();
    std::mt19937 rng(314159);
    std::uniform_real_distribution<double> coordinate(-20, 120), direction(-1, 1);
    std::vector<std::array<double, 3>> queries;
    for (int i = 0; i < 600; ++i)
        queries.push_back({{coordinate(rng), coordinate(rng), coordinate(rng)}});
    // Exact vertices and face centers exercise zero-distance and tied triangles.
    for (int i = 0; i < std::min(gen.triMesh.vNum, 60); ++i)
        queries.push_back({{gen.triMesh.v[i][0], gen.triMesh.v[i][1], gen.triMesh.v[i][2]}});
    struct Answer { int triangle; double distance; std::array<double, 3> point; };
    std::vector<Answer> reference;
    const auto start = std::chrono::steady_clock::now();
    for (auto& point : queries) {
        Answer best{-1, static_cast<double>(MAX_NUM2), {{0,0,0}}};
        for (int t = 0; t < gen.triMesh.eNum; ++t) {
            double closest[3];
            const double d = PointToTri(gen.triMesh.v[gen.triMesh.e[t][0]],
                gen.triMesh.v[gen.triMesh.e[t][1]], gen.triMesh.v[gen.triMesh.e[t][2]],
                point.data(), closest, best.distance);
            if (d < best.distance) {
                best.triangle = t; best.distance = d;
                std::copy(closest, closest + 3, best.point.begin());
            }
        }
        reference.push_back(best);
    }
    const auto middle = std::chrono::steady_clock::now();
    size_t evaluated = 0;
    for (size_t i = 0; i < queries.size(); ++i) {
        const auto best = gen.surfaceIndex.Nearest(queries[i].data(), [&](int t, double* closest) {
            ++evaluated;
            return PointToTri(gen.triMesh.v[gen.triMesh.e[t][0]], gen.triMesh.v[gen.triMesh.e[t][1]],
                gen.triMesh.v[gen.triMesh.e[t][2]], queries[i].data(), closest,
                std::numeric_limits<double>::infinity());
        });
        if (best.triangle != reference[i].triangle)
            std::cerr << "query " << queries[i][0] << ' ' << queries[i][1] << ' ' << queries[i][2]
                      << " brute " << reference[i].triangle << ':' << reference[i].distance
                      << " bvh " << best.triangle << ':' << best.distance << '\n';
        Require(best.triangle == reference[i].triangle, "Nearest triangle/tie mismatch at query " + std::to_string(i));
        Require(best.distance == reference[i].distance, "Nearest distance mismatch");
        Require(std::memcmp(best.point.data(), reference[i].point.data(), 3 * sizeof(double)) == 0,
                "Nearest point mismatch");
    }
    const auto end = std::chrono::steady_clock::now();
    std::cout << "triangles=" << gen.triMesh.eNum << " queries=" << queries.size()
              << " brute_seconds=" << std::chrono::duration<double>(middle-start).count()
              << " bvh_seconds=" << std::chrono::duration<double>(end-middle).count()
              << " triangle_evaluations=" << evaluated << "/" << queries.size()*gen.triMesh.eNum << '\n';

    std::vector<int> candidates;
    for (int i = 0; i < 600; ++i) {
        double point[3] = {coordinate(rng), coordinate(rng), coordinate(rng)};
        double dir[3] = {direction(rng), direction(rng), direction(rng)};
        if (i < 3) { dir[0] = dir[1] = dir[2] = 0; dir[i] = 1; }
        // Both unlimited lines (inside/outside) and finite segments (thickness).
        const double reach = i % 2 ? 16 / std::max(std::abs(dir[0]),
            std::max(std::abs(dir[1]), std::abs(dir[2]))) : std::numeric_limits<double>::infinity();
        gen.surfaceIndex.LineCandidates(point, dir, -reach, reach, candidates);
        std::vector<int> expected, actual;
        for (int t = 0; t < gen.triMesh.eNum; ++t) {
            double hit[3], alpha = 0;
            int code = Intersect(gen.triMesh.v[gen.triMesh.e[t][0]], gen.triMesh.v[gen.triMesh.e[t][1]],
                gen.triMesh.v[gen.triMesh.e[t][2]], point, dir, hit, alpha);
            if (code == 1 && alpha >= -reach && alpha <= reach) expected.push_back(t);
        }
        for (int t : candidates) {
            double hit[3], alpha = 0;
            int code = Intersect(gen.triMesh.v[gen.triMesh.e[t][0]], gen.triMesh.v[gen.triMesh.e[t][1]],
                gen.triMesh.v[gen.triMesh.e[t][2]], point, dir, hit, alpha);
            if (code == 1 && alpha >= -reach && alpha <= reach) actual.push_back(t);
        }
        Require(actual == expected, "BVH missed a line/segment intersection");
    }
}

static void Split(hexGen& gen, size_t index) {
    const int cell = gen.cutArray[index], level = gen.GetLevel(cell);
    gen.octreeArray[cell] = true;
    gen.cutArray.erase(gen.cutArray.begin() + index);
    for (int i = 0; i < 8; ++i) gen.cutArray.push_back(gen.Child(cell, level, i));
}

static void ValidatePartition(hexGen& gen, int resolutionDepth) {
    const int n = 1 << resolutionDepth;
    std::vector<int> levels(n*n*n, -1);
    Require(gen.leafNum == static_cast<int>(gen.cutArray.size()), "Leaf count mismatch");
    Require(std::unordered_set<int>(gen.cutArray.begin(), gen.cutArray.end()).size() == gen.cutArray.size(),
            "Duplicate leaf IDs");
    for (int cell : gen.cutArray) {
        const int level = gen.GetLevel(cell), size = n >> level;
        Require(size > 0, "Unexpected refinement beyond the test grid");
        int x,y,z; gen.OctreeidxToXYZ(cell,x,y,z,level);
        for (int i=x*size; i<(x+1)*size; ++i)
            for (int j=y*size; j<(y+1)*size; ++j)
                for (int k=z*size; k<(z+1)*size; ++k) {
                    int& occupied = levels[(i*n+j)*n+k];
                    Require(occupied == -1, "Overlapping leaves");
                    occupied = level;
                }
    }
    for (int value : levels) Require(value >= 0, "Hole in the leaf partition");
    for (int x=0; x<n; ++x) for (int y=0; y<n; ++y) for (int z=0; z<n; ++z)
        for (int dx=-1; dx<=1; ++dx) for (int dy=-1; dy<=1; ++dy) for (int dz=-1; dz<=1; ++dz) {
            const int xx=x+dx, yy=y+dy, zz=z+dz;
            if (xx<0 || yy<0 || zz<0 || xx>=n || yy>=n || zz>=n) continue;
            Require(std::abs(levels[(x*n+y)*n+z]-levels[(xx*n+yy)*n+zz]) <= 1,
                    "Touching leaves differ by more than one level");
        }
}

static void Balance() {
    for (int seed : {1, 7, 42}) {
        hexGen gen(10); // Sparse corners at the real 1025^3 lattice resolution.
        gen.octreeArray.resize(levelId[7], false);
        gen.cutArray.push_back(0);
        std::mt19937 random(seed);
        for (int i=0; i<45; ++i) {
            size_t index = random() % gen.cutArray.size();
            if (gen.GetLevel(gen.cutArray[index]) < 5) Split(gen,index);
        }
        const auto before = gen.cutArray.size();
        gen.StrongBalancedOctree();
        ValidatePartition(gen, 6);
        const auto balanced = gen.cutArray;
        gen.StrongBalancedOctree();
        Require(gen.cutArray == balanced, "Balancing is not idempotent");
        std::cout << "balance seed=" << seed << " before=" << before << " after=" << gen.cutArray.size() << '\n';
    }
    hexGen gen(10);
    for (int level=0; level<=10; ++level) {
        Require(gen.GetLevel(levelId[level]) == level, "Level lower bound mismatch");
        Require(gen.GetLevel(levelId[level+1]-1) == level, "Level upper bound mismatch");
    }
}

static void InteriorCube(hexGen& gen) {
    CubeSurface(gen);
    gen.octreeMesh.Initialize(32,64);
    gen.octreeMesh.vNum = 8; gen.octreeMesh.eNum = gen.leafNum = 1;
    gen.BOX_LENGTH_RATIO = 1;
    gen.START_POINT[0] = gen.START_POINT[1] = gen.START_POINT[2] = 0;
    for (int i=0; i<8; ++i) {
        gen.octreeMesh.e[0][i] = i;
        for (int a=0; a<3; ++a) gen.octreeMesh.v[i][a] = 40 + gen.triMesh.v[i][a] * 0.2;
    }
}

static void Projection() {
    ProjectionOptions options;
    options.maxIterations = 5; options.stagnationChecks = 2;
    ProjectionControl control(options);
    Require(control.Check(1,2,0.1,10) == ProjectionStatus::Running, "Premature stop");
    Require(control.Check(2,2,0.1,10) == ProjectionStatus::Running, "Premature stagnation");
    Require(control.Check(3,2,0.1,10) == ProjectionStatus::Stagnated, "Missing stagnation stop");
    Require(control.Check(5,2,0.1,10) == ProjectionStatus::IterationLimit, "Missing iteration limit");
    Require(control.Check(5,0,0.6,0) == ProjectionStatus::Converged, "Convergence at budget boundary");
    Require(control.Check(1,0,std::numeric_limits<double>::quiet_NaN(),0) == ProjectionStatus::InvalidGeometry,
            "Nonfinite quality accepted");
    options.checkEvery = 0;
    bool rejected = false;
    try { ProjectionControl invalid(options); } catch (const std::invalid_argument&) { rejected = true; }
    Require(rejected, "Invalid stopping settings accepted");

    hexGen gen(5);
    InteriorCube(gen);
    options = ProjectionOptions();
    options.maxIterations = 1;
    std::srand(1);
    auto result = gen.ProjectToIsoSurface("limit-test.vtk",options);
    Require(result.status == ProjectionStatus::IterationLimit && result.iterations == 1,
            "Projection did not stop at a non-checkpoint iteration limit");
    std::remove("limit-test.vtk");

    hexGen convergence(5);
    InteriorCube(convergence);
    // Exercise successful stopping at an actual measured checkpoint with a
    // deliberately loose surface tolerance, independent of optimizer speed.
    options.targetScaledJacobian = 0.02;
    options.surfaceTolerance = 40;
    options.maxIterations = 1;
    options.checkEvery = 100;
    options.stagnationChecks = 100;
    std::srand(1);
    result = convergence.ProjectToIsoSurface("convergence-test.vtk",options);
    Require(result.status == ProjectionStatus::Converged, "Cube projection failed to converge");
    Require(result.minScaledJacobian > options.targetScaledJacobian &&
            result.maxSurfaceDistance < options.surfaceTolerance, "Invalid convergence metrics");
    std::remove("convergence-test.vtk");
    std::remove("finalMesh.vtk");
}

static void Pipeline(const char* input) {
    if (!input) {
        hexGen surface(3);
        CubeSurface(surface);
        FILE* file = std::fopen("pipeline-input.raw", "w");
        Require(file != nullptr, "Cannot create pipeline input");
        std::fprintf(file, "8 12\n");
        for (int i=0; i<8; ++i)
            std::fprintf(file, "%lf %lf %lf\n", surface.triMesh.v[i][0], surface.triMesh.v[i][1], surface.triMesh.v[i][2]);
        for (int i=0; i<12; ++i)
            std::fprintf(file, "%d %d %d\n", surface.triMesh.e[i][0], surface.triMesh.e[i][1], surface.triMesh.e[i][2]);
        std::fclose(file);
        input = "pipeline-input.raw";
    }
    hexGen gen(4);
    gen.InitializeOctree(input, "pipeline-surface.vtk");
    gen.ConstructOctree();
    ValidatePartition(gen, 4);
    gen.OutputOctree("pipeline-octree.vtk");
    gen.DualFullHexMeshExtraction("pipeline-full.vtk");
    std::srand(1);
    gen.RemoveOutsideElement("pipeline-interior.vtk");
    Require(gen.leafNum > 0, "Pipeline produced no interior mesh");
    ProjectionOptions options;
    options.maxIterations = 2;
    const auto result = gen.ProjectToIsoSurface("pipeline-projected.vtk", options);
    Require(result.status == ProjectionStatus::IterationLimit || result.status == ProjectionStatus::Converged,
            "Pipeline projection failed");
    for (int i=0; i<gen.octreeMesh.vNum; ++i)
        for (int axis=0; axis<3; ++axis)
            Require(std::isfinite(gen.octreeMesh.v[i][axis]), "Nonfinite pipeline vertex");
    for (int i=0; i<gen.octreeMesh.eNum; ++i)
        for (int j=0; j<8; ++j)
            Require(gen.octreeMesh.e[i][j] >= 0 && gen.octreeMesh.e[i][j] < gen.octreeMesh.vNum,
                    "Invalid output connectivity");
    std::cout << "pipeline vertices=" << gen.octreeMesh.vNum << " elements=" << gen.octreeMesh.eNum << '\n';
    for (const char* name : {"pipeline-input.raw", "pipeline-surface.vtk", "pipeline-octree.vtk",
        "pipeline-full.vtk", "pipeline-interior.vtk", "pipeline-projected.vtk"}) std::remove(name);
}

int main(int argc, char** argv) try {
    Require(argc >= 2, "Missing test mode");
    std::string mode = argv[1];
    if (mode == "geometry") {
        hexGen cube(5); CubeSurface(cube); Geometry(cube);
        for (int i=2; i<argc; ++i) {
            hexGen gen(10);
            gen.ReadRawData(argv[i], "regression-surface.vtk");
            Geometry(gen);
        }
        std::remove("regression-surface.vtk");
    } else if (mode == "balance") Balance();
    else if (mode == "projection") Projection();
    else if (mode == "pipeline") Pipeline(argc > 2 ? argv[2] : nullptr);
    else throw std::runtime_error("Unknown test mode");
    std::cout << mode << " passed\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
