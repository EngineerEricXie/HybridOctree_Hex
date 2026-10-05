[![DOI](https://img.shields.io/badge/DOI-10.1016/j.jocs.2024.102278-blue)](https://doi.org/10.1016/j.jocs.2024.102278)
# HybridOctree_Hex
Please download the latest release of the software.

HybridOctree_Hex received the best technical poster award in [the 2024 International Meshing Roundtable](https://internationalmeshingroundtable.com/awards/).

Please also feel free to check out other works on grid-based hexahedral mesh generation:

1. https://github.com/CMU-CBML/Element-Saving-Hexahedral-3-Refinement-Templates

## Build

Use an out-of-source build with CMake 3.16 or later and a C++11 compiler:

```sh
cmake -S HybridOctree_Hex -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release --config Release --parallel
```

Single-configuration generators now default to `Release` when no build type is
specified. Explicit `Debug` and other configurations are preserved. The checked-in
generated Makefiles/cache refer to a Windows installation; generate a fresh build
directory on other machines. `Main.cpp` reads `model.raw` from the working directory.

## Bounded projection, sparse balancing and BVH queries

The three follow-up optimizations are implemented:

- `ProjectToIsoSurface` returns a `ProjectionResult` with a stopping reason,
  iteration count, minimum scaled Jacobian and maximum surface distance. The
  quality threshold stops increasing at the requested target. Iteration counting
  is independent of the smoothing loop and includes the last partial checkpoint.
  `Main.cpp` returns 0 on success, 2 when projection stops without meeting the
  target, and 1 on an exception. The current finite checkpoint is written to
  `projHex.vtk`; `finalMesh.vtk` is refreshed only after convergence.
- `StrongBalancedOctree` stores occupied integer lattice vertices and performs
  iterative passes. It finds touching coarse leaves even when a fine corner lies
  inside a coarse face or edge, including domain boundaries. Siblings are refined
  together for the extraction templates. The fixed `(1025)^3` vector allocation
  (about 24.07 GiB of headers per old invocation) and recursion are removed.
  Leaf lists now contain only active IDs, with no 100-million-entry padding.
  Levels are computed from ID ranges, eliminating the dense level table and
  correctly handling the deepest level.
- A reusable `TriangleBVH` filters thickness segments and inside/outside lines,
  and performs nearest-surface queries for removal and projection. Existing
  triangle predicates remain the narrow phase. Candidate IDs are sorted, and
  nearest-distance ties choose the lowest input triangle ID. The BVH is built
  lazily and invalidated when the surface is read again.

Projection defaults are in `ProjectionControl.h`; configure `projectionOptions`
in `Main.cpp` or pass options directly to `ProjectToIsoSurface`:

| Option | Default | Meaning |
| --- | ---: | --- |
| `targetScaledJacobian` | 0.53 | Converged mesh must have minimum SJ strictly above this value. |
| `surfaceTolerance` | 1e-6 | Maximum surface distance must be below this value, in normalized [0, 100] coordinates. |
| `maxIterations` | 100000 | Hard limit on gradient iterations. |
| `checkEvery` | 1000 | Quality, distance, smoothing and checkpoint interval; the last iteration is always checked. |
| `stagnationChecks` | 20 | Stop after this many consecutive checkpoints without significant improvement. |
| `progressTolerance` | 1e-8 | Absolute improvement threshold for SJ and distance; a reduction in bad-element count also resets stagnation. |

Nonfinite updated coordinates/quality stop with `InvalidGeometry`. A finite
iteration budget does not guarantee that the optimizer can attain the requested
quality on every input. Inspect the returned status before using an output as a
converged result.

### Validation and measurements

```sh
cmake -S HybridOctree_Hex -B build/release -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build/release --parallel
ctest --test-dir build/release --output-on-failure
# Print the brute-force/BVH query timings and triangle-evaluation counts:
ctest --test-dir build/release -R geometry_queries -V
```

Five CTest cases pass in Release and with AddressSanitizer/UndefinedBehaviorSanitizer:
geometry queries, sparse balancing, projection stopping, and full-stage cube/bone
smoke tests. Geometry checks compare triangle IDs, distances and closest-point
bytes against the full scan on cube, bone and sphereInShell, plus full-line and
finite-segment intersection sets. Balance tests verify complete coverage, no
overlap/duplicate leaves, all 26 kinds of voxel adjacency satisfying 2:1 balance,
and idempotence on three adaptive partitions at a depth-10 lattice resolution.
Projection tests cover successful stopping (using deliberately loose fixture
targets), stagnation, invalid options, nonfinite quality and a limit that is not a
regular checkpoint. Full-stage tests use depth 4 and two projection iterations;
they validate finite coordinates and connectivity, not default-target convergence
of all production models. Leak detection is disabled for the sanitizer run because
the pre-existing mesh ownership/destructor leaks remain outside this change.

On the same i9-14900KF/GCC 11.4 machine, bone initialization now takes about
**0.047 s and 151 MiB peak RSS**, versus **0.661 s and 2.04 GiB** for the previous
byte-level-table version (three-run medians, identical `-O3` flags). The original
pre-optimization initialization required about 5.46 GiB. A dense refinement bit
array still occupies about 146 MiB at depth 10, and construction releases it after
balancing.

For 660 nearest-surface queries per model, the BVH evaluates **31,521 instead of
7,978,080 triangles on bone**, and **14,132 instead of 4,051,080 on sphereInShell**.
Measured nearest-query speedups are approximately 50x and 30x respectively;
small meshes such as the 12-triangle cube may be slower due to traversal overhead.
These are stage measurements, not whole-application speedup claims.

### Numerical behavior changes

Unqualified `abs(double)` resolved to the integer overload with GCC, truncating
distances and direction components. Calls now use `std::abs`, retaining fractions;
this is necessary for valid distance bounds during BVH pruning. Broad-phase
filtering also avoids retrying rays because of unrelated parallel triangles.
Ambiguous candidate intersections still retry, up to 128 attempts before an
explicit error. Together with corrected deepest-level IDs and coarse/fine
balancing, these fixes mean full meshes need not be byte-identical to the old
implementation. The BVH differential tests use the corrected floating-point
predicates on both sides. Curvature preprocessing remains byte-identical.

Remaining performance work includes vertex deduplication, face/element adjacency
indexes, and hierarchical cell construction. The unused level-9 scan is skipped
(19,173,961 rather than 153,391,689 cell visits at default depth), but active-level
cell/triangle scans remain. The gradient optimizer itself was not redesigned.

## Initial preprocessing measurements (2026-10-05)

The following records the first optimization pass, before the changes above.
Its curvature calculation and traversal order were preserved:

- `ReadRawData` builds a vertex-to-triangle index and checks only incident
  triangles when accumulating curvature. The previous implementation searched
  every remaining triangle for each vertex of each triangle. For ordinary meshes
  with bounded vertex valence, candidate work is approximately linear instead of
  quadratic; unusually high valence can still cause quadratic work. The index
  requires O(vertices + triangles) additional temporary storage.
- The intermediate version stored `getLevel` as `unsigned char` instead of `int`.
  At the default `VOXEL_SIZE=10`, this removes approximately 3.43 GiB of allocation
  and zero-initialization on platforms with four-byte integers. Table size and
  values were preserved, including the original zero-filled deepest level. The
  current implementation replaces this table with `GetLevel` as described above.

Measured on an Intel Core i9-14900KF under WSL/Linux, GCC 11.4, `-O3 -DNDEBUG`,
against commit `dfcb3db`. Both versions used the same compiler flags. Values are
the median of three separate process runs, alternating baseline/candidate order.
These are **surface reading, normalization, curvature calculation and VTK output**
times; they exclude octree allocation/construction and subsequent mesh stages.

| Input | Triangles | Before | After | Stage speedup |
| --- | ---: | ---: | ---: | ---: |
| sphereInShell | 6,138 | 0.0346 s | 0.00410 s | 8.4x |
| bone | 12,088 | 0.1268 s | 0.00838 s | 15.1x |
| bunny | 22,490 | 0.4516 s | 0.01719 s | 26.3x |
| bumpy_torus | 30,558 | 0.8291 s | 0.02441 s | 34.0x |
| fertility | 52,456 | 2.4892 s | 0.03762 s | 66.2x |

Including default octree allocation, `InitializeOctree` on bone decreased from
1.953 s to 0.589 s (3.3x). Peak process RSS decreased from approximately 5.46 GiB
to 2.04 GiB (63%). These results are environment dependent and are not a claim
about total mesh-generation speed.

Validation compared VTK bytes and the unrounded binary `double` curvature arrays:
all five meshes and six synthetic cases matched exactly. Synthetic cases cover
boundary edges, vertex-only adjacency, nonmanifold edges, duplicate faces,
repeated vertex IDs, and mixed triangle order. The initialization benchmark also
checks old table values or current level-range boundaries outside the timed section. A Release build passed;
the source's missing `<cmath>` dependency was added to make it compile with GCC.

Reproduce the measurements and equivalence checks with Python 3, GCC and GNU time:

```sh
baseline_dir="$(mktemp -d)"
git archive dfcb3db HybridOctree_Hex | tar -x -C "$baseline_dir"
python3 tools/benchmark_curvature.py "$baseline_dir/HybridOctree_Hex" \
  "input boundaries/sphereInShell_tri.raw" "input boundaries/bone_tri.raw" \
  "input boundaries/bunny_tri.raw" "input boundaries/bumpy_torus_tri.raw" \
  "input boundaries/fertility_tri.raw"
python3 tools/benchmark_curvature.py "$baseline_dir/HybridOctree_Hex" \
  "input boundaries/bone_tri.raw" --mode initialize
```

The benchmark uses a temporary driver to access preprocessing directly; it does
not change the application API. Initialization mode needs enough RAM for the
baseline's approximately 5.5 GiB allocation. The synthetic cases test equivalence
with existing behavior, not the geometric validity of degenerate input.

# Citation
```angular2html
@article{tong2024hybridoctree_hex,
  title={HybridOctree\_Hex: Hybrid octree-based adaptive all-hexahedral mesh generation with Jacobian control},
  author={Tong, Hua and Halilaj, Eni and Zhang, Yongjie Jessica},
  journal={Journal of Computational Science},
  volume={78},
  pages={102278},
  year={2024},
  publisher={Elsevier}
}
```

# Mesh Statistics:
|Model|#Vert|#Elem|Min SJ|
| :--: | :--: | :--: | :--: |
|airplane1|29899|24318|0.57|
|airplane2|29241|23646|0.61|
|ant2|57833|48120|0.60|
|armadillo|326447|282201|0.54|
|armchair|75577|65280|0.60|
|bimba|370664|330619|0.57|
|bird|25785|21067|0.58|
|blade|103933|88351|0.56|
|bone|10356|8619|0.61|
|botijo|161860|139327|0.57|
|bottle1|36091|30145|0.56|
|bottle2|89357|78366|0.58|
|bumpy_sphere|25466|21646|0.58|
|bumpy_torus|249852|219832|0.56|
|bunny|26375|21695|0.57|
|buste|71068|61441|0.56|
|camel|124573|105417|0.57|
|camille_hand|114453|100023|0.59|
|carter|190871|160122|0.55|
|chair|232460|192814|0.57|
|chair1|50650|41034|0.58|
|chinese_lion|351402|309525|0.55|
|cup|386677|342009|0.58|
|Cup1|89032|73266|0.58|
|dancer|58353|47994|0.56|
|dancer2|43774|35142|0.57|
|dancing_children|457139|396964|0.54|
|david|319465|282957|0.56|
|deformed_armadillo|43216|34939|0.56|
|dente|129398|114812|0.58|
|dilo|94803|79929|0.60|
|dino|113019|94303|0.55|
|dino2|75158|62562|0.56|
|dragonstand2|62576|50853|0.56|
|dragon_stand|240896|206695|0.55|
|dtorus|19515|15872|0.58|
|duck|65685|58538|0.58|
|duck2|218549|193450|0.58|
|eight|30157|25264|0.61|
|elephant|59595|49242|0.57|
|elk|355959|313725|0.56|
|eros|75907|65798|0.57|
|fertility|45514|38475|0.58|
|fish1|64760|53611|0.58|
|fish2|51367|42853|0.60|
|foot|64617|56468|0.58|
|gargoyle|273704|236689|0.55|
|genus3|54592|45176|0.57|
|glass1|42870|34034|0.57|
|glass2|14791|11602|0.62|
|grayloc|1102475|986401|0.56|
|greek_sculpture|187058|162674|0.57|
|hand|88887|75770|0.57|
|head|62782|55038|0.55|
|head1|153260|134611|0.56|
|head2|361054|327103|0.57|
|holes3|41611|36132|0.61|
|homer|166322|146056|0.57|
|horse|102832|87595|0.57|
|human1|148026|128154|0.60|
|human2|141982|119792|0.54|
|human3|187747|163852|0.55|
|igea|74420|64303|0.56|
|insect|108197|89154|0.54|
|isidore_horse|209974|182124|0.54|
|Kiss|116510|100375|0.55|
|kitten|38138|32258|0.55|
|lion_recon|112847|95251|0.57|
|master_cylinder|286569|247194|0.53|
|max|55222|47687|0.56|
|moai|145934|127950|0.54|
|mouse|99657|85890|0.53|
|oil_pump|233702|196455|0.54|
|oni|135607|117013|0.55|
|pear|19308|16122|0.58|
|pensatore|380202|340985|0.56|
|pierrot|170778|151683|0.55|
|pig|84664|72571|0.56|
|ramses|44790|37993|0.59|
|red_circular_box|351881|313866|0.56|
|retinal|19021|15950|0.58|
|rocker|93607|80266|0.54|
|rolling_stage|343540|297087|0.53|
|santa|93459|78597|0.54|
|screwdriver|147661|126098|0.53|
|sediapatch|56036|47445|0.55|
|sphereInShell|128626|108381|0.62|
|sphinx|206944|183216|0.55|
|teaport|20093|16726|0.57|
|thai_statue|64764|53831|0.55|
|toy1|20274|16660|0.58|
|toy2|15818|13268|0.58|
|uu-memento|34335|27228|0.58|
|venus|25794|21895|0.58|
|woodenfish|254517|214410|0.55|
|wrench|56712|45382|0.54|
