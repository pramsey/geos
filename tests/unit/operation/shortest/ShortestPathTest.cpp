//
// Test Suite for geos::operation::shortest::ShortestPath class.

// tut
#include <tut/tut.hpp>
// geos
#include <geos/operation/shortest/ShortestPath.h>
#include <geos/io/WKTReader.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Curve.h>
#include <geos/geom/LineString.h>
#include <geos/geom/GeometryCollection.h>
// std
#include <memory>
#include <string>
#include <vector>

namespace tut {
//
// Test Group
//

struct test_shortest_data {
    geos::io::WKTReader wktreader;

    test_shortest_data()
        : wktreader()
    {
    }

    std::unique_ptr<geos::geom::Geometry>
    readWKT(const std::string& inputWKT)
    {
        return std::unique_ptr<geos::geom::Geometry>(wktreader.read(inputWKT));
    }

    std::vector<const geos::geom::Curve*>
    toCurves(const geos::geom::Geometry* geom) {
        std::vector<const geos::geom::Curve*> curves;
        if (const geos::geom::GeometryCollection* gc = dynamic_cast<const geos::geom::GeometryCollection*>(geom)) {
            for (std::size_t i = 0; i < gc->getNumGeometries(); ++i) {
                if (const geos::geom::Curve* c = dynamic_cast<const geos::geom::Curve*>(gc->getGeometryN(i))) {
                    curves.push_back(c);
                }
            }
        } else if (const geos::geom::Curve* c = dynamic_cast<const geos::geom::Curve*>(geom)) {
            curves.push_back(c);
        }
        return curves;
    }

    geos::geom::CoordinateXY pt(double x, double y) {
        return geos::geom::CoordinateXY(x, y);
    }
};

typedef test_group<test_shortest_data> group;
typedef group::object object;

group test_shortest_group("geos::operation::shortest::ShortestPath");

//
// Test Cases
//

// Test 1: single edge, direct path
template<> template<>
void object::test<1>()
{
    auto geom = readWKT("LINESTRING(0 0, 10 0)");
    auto curves = toCurves(geom.get());

    std::vector<std::size_t> result;
    geos::operation::shortest::ShortestPath::shortestPath(
        curves, pt(0,0), pt(10,0), result);

    ensure_equals(result.size(), 1u);
    ensure_equals(result[0], 1u);
}

// Test 2: two-edge chain — both edges in path with correct ordering
template<> template<>
void object::test<2>()
{
    auto geom = readWKT("MULTILINESTRING((0 0, 10 0), (10 0, 20 0))");
    auto curves = toCurves(geom.get());

    std::vector<std::size_t> result;
    geos::operation::shortest::ShortestPath::shortestPath(
        curves, pt(0,0), pt(20,0), result);

    ensure_equals(result.size(), 2u);
    ensure_equals(result[0], 1u);
    ensure_equals(result[1], 2u);
}

// Test 3: triangle — direct edge is shorter than two-edge route
// Edges:
//   0: (0,0)→(10,0)  length 10
//   1: (10,0)→(5,10) length ~11.18
//   2: (5,10)→(0,0)  length ~11.18
// Shortest path from (0,0) to (5,10): direct via edge 2 (~11.18)
// versus via edge 0+1 (~21.18). Should pick edge 2 (position 1).
template<> template<>
void object::test<3>()
{
    auto geom = readWKT("MULTILINESTRING((0 0, 10 0), (10 0, 5 10), (5 10, 0 0))");
    auto curves = toCurves(geom.get());

    std::vector<std::size_t> result;
    geos::operation::shortest::ShortestPath::shortestPath(
        curves, pt(0,0), pt(5,10), result);

    ensure_equals(result.size(), 3u);

    // Only one edge should be selected
    int count = 0;
    for (std::size_t r : result) if (r > 0) count++;
    ensure_equals(count, 1);

    // Edge 2 (direct (0,0)↔(5,10)) is shortest
    ensure_equals(result[2], 1u);
    ensure_equals(result[0], 0u);
    ensure_equals(result[1], 0u);
}

// Test 4: snap to nearest node — start/end not exactly on node
template<> template<>
void object::test<4>()
{
    auto geom = readWKT("MULTILINESTRING((0 0, 10 0), (10 0, 20 0))");
    auto curves = toCurves(geom.get());

    std::vector<std::size_t> result;
    // Start near (0,0), end near (20,0)
    geos::operation::shortest::ShortestPath::shortestPath(
        curves, pt(0.1, 0.1), pt(19.9, 0.0), result);

    ensure_equals(result.size(), 2u);
    ensure_equals(result[0], 1u);
    ensure_equals(result[1], 2u);
}

// Test 5: start and end snap to same node — no path
template<> template<>
void object::test<5>()
{
    auto geom = readWKT("MULTILINESTRING((0 0, 10 0), (10 0, 20 0))");
    auto curves = toCurves(geom.get());

    std::vector<std::size_t> result;
    geos::operation::shortest::ShortestPath::shortestPath(
        curves, pt(0,0), pt(0,0), result);

    ensure_equals(result.size(), 2u);
    ensure_equals(result[0], 0u);
    ensure_equals(result[1], 0u);
}

// Test 6: empty input
template<> template<>
void object::test<6>()
{
    std::vector<const geos::geom::Curve*> curves;
    std::vector<std::size_t> result;
    geos::operation::shortest::ShortestPath::shortestPath(
        curves, pt(0,0), pt(10,0), result);
    ensure(result.empty());
}

// Test 7: null and empty geometry entries
template<> template<>
void object::test<7>()
{
    auto factory = geos::geom::GeometryFactory::create();
    auto empty = factory->createLineString();

    std::vector<const geos::geom::Curve*> curves;
    curves.push_back(nullptr);
    curves.push_back(empty.get());

    std::vector<std::size_t> result;
    geos::operation::shortest::ShortestPath::shortestPath(
        curves, pt(0,0), pt(10,0), result);

    ensure_equals(result.size(), 2u);
    ensure_equals(result[0], 0u);
    ensure_equals(result[1], 0u);
}

// Test 8: disconnected graph — no path exists
template<> template<>
void object::test<8>()
{
    auto geom = readWKT("MULTILINESTRING((0 0, 10 0), (20 0, 30 0))");
    auto curves = toCurves(geom.get());

    std::vector<std::size_t> result;
    geos::operation::shortest::ShortestPath::shortestPath(
        curves, pt(0,0), pt(30,0), result);

    ensure_equals(result.size(), 2u);
    ensure_equals(result[0], 0u);
    ensure_equals(result[1], 0u);
}

// Test 9: longestShortestPath on a 3-edge chain — all 3 edges on path
template<> template<>
void object::test<9>()
{
    auto geom = readWKT("MULTILINESTRING((0 0, 10 0), (10 0, 20 0), (20 0, 30 0))");
    auto curves = toCurves(geom.get());

    std::vector<std::size_t> result;
    geos::operation::shortest::ShortestPath::longestShortestPath(curves, result);

    ensure_equals(result.size(), 3u);

    int count = 0;
    for (std::size_t r : result) if (r > 0) count++;
    ensure_equals(count, 3);

    // Positions should be 1,2,3 in some order (all edges on path)
    for (std::size_t r : result) {
        ensure(r >= 1u && r <= 3u);
    }
}

// Test 10: longestShortestPath on a star graph
// 4 arms from (0,0): lengths 10, 10, 10, 10
// Longest path = 2 edges (one arm through center to opposite arm)
template<> template<>
void object::test<10>()
{
    auto geom = readWKT("MULTILINESTRING((0 0, 10 0), (0 0, 0 10), (0 0, -10 0), (0 0, 0 -10))");
    auto curves = toCurves(geom.get());

    std::vector<std::size_t> result;
    geos::operation::shortest::ShortestPath::longestShortestPath(curves, result);

    ensure_equals(result.size(), 4u);

    int count = 0;
    for (std::size_t r : result) if (r > 0) count++;
    ensure_equals(count, 2);
}

// Test 11: longestShortestPath on empty input
template<> template<>
void object::test<11>()
{
    std::vector<const geos::geom::Curve*> curves;
    std::vector<std::size_t> result;
    geos::operation::shortest::ShortestPath::longestShortestPath(curves, result);
    ensure(result.empty());
}

// Test 12: shortest path reverse direction — same result regardless of start/end order
template<> template<>
void object::test<12>()
{
    auto geom = readWKT("MULTILINESTRING((0 0, 10 0), (10 0, 20 0))");
    auto curves = toCurves(geom.get());

    std::vector<std::size_t> resultFwd, resultRev;
    geos::operation::shortest::ShortestPath::shortestPath(
        curves, pt(0,0), pt(20,0), resultFwd);
    geos::operation::shortest::ShortestPath::shortestPath(
        curves, pt(20,0), pt(0,0), resultRev);

    // Both directions should select the same 2 edges (though positions may differ)
    int countFwd = 0, countRev = 0;
    for (std::size_t r : resultFwd) if (r > 0) countFwd++;
    for (std::size_t r : resultRev) if (r > 0) countRev++;
    ensure_equals(countFwd, 2);
    ensure_equals(countRev, 2);
}

} // namespace tut
