This is a C/C++ computational geometry library (GEOS 3.15-dev), a C++ port of the Java
[JTS Topology Suite](https://github.com/locationtech/jts/) v1.18.0.

## Building

Use `_build` as the build directory. From `_build`: `cmake .. && make`.
`CMAKE_BUILD_TYPE` options: `Release`, `Debug`, `ASAN`, `UBSAN`.

Build outputs:

* `./bin/test_geos_unit` — unit test runner; `-h` for help
  * `./bin/test_geos_unit geos::triangulate::Delaunay` — single suite
  * `./bin/test_geos_unit geos::triangulate::Delaunay 4` — single test
* `./bin/test_xmltester` — XML test runner (`-v` summary, `-v -v` with geometry dumps)
* `./bin/geosop` — CLI harness; `-h` lists all available operations; useful for
  manually verifying that an operation produces the expected output

## Source Layout

Headers in `include/geos/`, sources in `src/`, both sharing the same module hierarchy:

| Module | Purpose |
|---|---|
| `geom` | Core geometry types (Point, LineString, Polygon, …) |
| `algorithm` | Core geometric algorithms (area, distance, orientation, …) |
| `operation` | High-level operations (buffer, overlay, distance, cluster, …) |
| `index` | Spatial indexes (STRtree, KdTree, quadtree, …) |
| `io` | WKT, WKB, GeoJSON readers/writers |
| `noding` | Noding and segment-string operations |
| `triangulate` | Delaunay triangulation and Voronoi diagrams |
| `coverage` | Coverage validation and cleaning |
| `precision` | Coordinate precision reduction |
| `simplify` | Douglas-Peucker and related simplification |
| `linearref` | Linear referencing along line geometries |
| `dissolve` | Line and geometry dissolution |
| `edgegraph` | Half-edge graph data structure |
| `geomgraph` | Geometric graph (used internally by overlay) |
| `planargraph` | Planar graph algorithms |
| `shape` | Geometric shape factories and generators |
| `math` | Double-double precision arithmetic |
| `util` | Assert, Interrupt, Profiler utilities |

The C API lives in `capi/` and has a stable ABI across releases. The C++ API has no
ABI stability guarantees and the library is intentionally renamed on each minor release.
Do not remove or change the signature of any existing C API function.

## JTS Relationship

GEOS is a port of JTS. **Do not rename classes, methods, variables, or members** — it
makes porting JTS updates harder. New algorithms should be prototyped in JTS first and
then ported to GEOS. Do not file JTS tickets; that is a human decision.

When a commit ports a JTS change, include the source JTS PR or commit URL in the
commit message body (see Commit Messages below).

## C++ Standard

C++17 is the minimum. C++17 features (`std::optional`, structured bindings,
`if constexpr`, etc.) are fine to use.

## Code Style

Full guidelines are in `DEVELOPER-NOTES.md`. Key points for quick reference:

* Run `tools/astyle.sh <file>` before committing (Stroustrup style, 4-space indent,
  120-char max line).
* Use `#pragma once` in all headers.
* Add `GEOS_DLL` to every public-API class declaration.
* Prefer `const Geometry&` for read-only geometry arguments; use `Geometry*` only
  when null is a valid input.
* Raw pointers and references are non-owning by convention; use `std::unique_ptr<T>`
  for ownership.

## Writing Tests

Unit tests use the TUT framework in `tests/unit/`, mirroring the `src/` directory
structure. New test files must be named `*Test.cpp` to be auto-discovered.

```cpp
// Test Suite for geos::foo::Bar
#include <tut/tut.hpp>
#include <geos/foo/Bar.h>
#include <geos/io/WKTReader.h>

namespace tut {
struct test_bar_data {
    geos::geom::GeometryFactory::Ptr factory_;
    geos::io::WKTReader reader_;
    test_bar_data() : factory_(GeometryFactory::create()), reader_(factory_.get()) {}
};
typedef test_group<test_bar_data> group;
typedef group::object object;
group test_bar_group("geos::foo::Bar");  // suite name used on the test_geos_unit CLI

template<> template<>
void object::test<1>() { ... }
```

## Commit Messages

* Reference the GitHub issue: `Fix crash in overlay (GH-1234)`.
* For JTS ports, add the JTS PR or commit URL on a second line:

```
Port minimum spanning tree from JTS (GH-1099)
https://github.com/locationtech/jts/pull/999
```

## NEWS.md

Add an entry under the current version for: net-new features, behaviour changes,
JTS ports with observable effect, and (on stable branches) every bug fix.
Format: `- Description (GH-NNNN, Author Name)`.
No entry needed for refactors, test additions, or build-system changes with no
user-visible effect.

## Active Development

Curved geometry support (`CircularString`, `CompoundCurve`, `CurvePolygon`, etc.)
is under active development in 3.15. When touching any code path that handles
geometry types, check whether nearby code already dispatches on curve types and
handle or explicitly exclude them as appropriate.
