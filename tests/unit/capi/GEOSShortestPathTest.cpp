//
// Test Suite for C-API GEOSShortestPath and GEOSLongestShortestPath

#include <tut/tut.hpp>
// geos
#include <geos_c.h>

#include <vector>

#include "capi_test_utils.h"

namespace tut {
//
// Test Group
//

struct test_capigeosshortestpath_data : public capitest::utility {
};

typedef test_group<test_capigeosshortestpath_data> group;
typedef group::object object;

group test_capigeosshortestpath_group("capi::GEOSShortestPath");

//
// Test Cases
//

// Simple two-edge chain: full path selected
template<>
template<>
void object::test<1> ()
{
    constexpr int size = 2;
    GEOSGeometry* geoms[size];
    geoms[0] = GEOSGeomFromWKT("LINESTRING(0 0, 10 0)");
    geoms[1] = GEOSGeomFromWKT("LINESTRING(10 0, 20 0)");

    size_t* result = GEOSShortestPath(geoms, size, 0, 0, 20, 0);

    ensure(nullptr != result);
    ensure_equals(result[0], 1u);
    ensure_equals(result[1], 2u);

    GEOSFree(result);
    for (auto& g : geoms) GEOSGeom_destroy(g);
}

// Triangle: direct edge is shorter
template<>
template<>
void object::test<2> ()
{
    constexpr int size = 3;
    GEOSGeometry* geoms[size];
    geoms[0] = GEOSGeomFromWKT("LINESTRING(0 0, 10 0)");   // len 10
    geoms[1] = GEOSGeomFromWKT("LINESTRING(10 0, 5 10)");  // len ~11.18
    geoms[2] = GEOSGeomFromWKT("LINESTRING(5 10, 0 0)");   // len ~11.18

    // From (0,0) to (5,10): direct edge (index 2) is shortest
    size_t* result = GEOSShortestPath(geoms, size, 0, 0, 5, 10);

    ensure(nullptr != result);

    int count = 0;
    for (int i = 0; i < size; ++i)
        if (result[i] > 0) count++;
    ensure_equals(count, 1);
    ensure_equals(result[2], 1u);
    ensure_equals(result[0], 0u);
    ensure_equals(result[1], 0u);

    GEOSFree(result);
    for (auto& g : geoms) GEOSGeom_destroy(g);
}

// Disconnected graph: no path
template<>
template<>
void object::test<3> ()
{
    constexpr int size = 2;
    GEOSGeometry* geoms[size];
    geoms[0] = GEOSGeomFromWKT("LINESTRING(0 0, 10 0)");
    geoms[1] = GEOSGeomFromWKT("LINESTRING(20 0, 30 0)");

    size_t* result = GEOSShortestPath(geoms, size, 0, 0, 30, 0);

    ensure(nullptr != result);
    ensure_equals(result[0], 0u);
    ensure_equals(result[1], 0u);

    GEOSFree(result);
    for (auto& g : geoms) GEOSGeom_destroy(g);
}

// Non-curve geometries are ignored
template<>
template<>
void object::test<4> ()
{
    constexpr int size = 3;
    GEOSGeometry* geoms[size];
    geoms[0] = GEOSGeomFromWKT("LINESTRING(0 0, 10 0)");
    geoms[1] = GEOSGeomFromWKT("POINT(5 5)");
    geoms[2] = GEOSGeomFromWKT("LINESTRING(10 0, 20 0)");

    size_t* result = GEOSShortestPath(geoms, size, 0, 0, 20, 0);

    ensure(nullptr != result);
    ensure_equals(result[0], 1u);
    ensure_equals(result[1], 0u); // POINT ignored
    ensure_equals(result[2], 2u);

    GEOSFree(result);
    for (auto& g : geoms) GEOSGeom_destroy(g);
}

// Empty input returns NULL
template<>
template<>
void object::test<5> ()
{
    size_t* result = GEOSShortestPath(nullptr, 0, 0, 0, 10, 0);
    ensure(nullptr == result);
}

// Null array entries handled
template<>
template<>
void object::test<6> ()
{
    constexpr int size = 2;
    GEOSGeometry* geoms[size];
    geoms[0] = GEOSGeomFromWKT("LINESTRING(0 0, 10 0)");
    geoms[1] = nullptr;

    size_t* result = GEOSShortestPath(geoms, size, 0, 0, 10, 0);

    ensure(nullptr != result);
    ensure_equals(result[0], 1u);
    ensure_equals(result[1], 0u);

    GEOSFree(result);
    GEOSGeom_destroy(geoms[0]);
}

// LongestShortestPath: 3-edge chain selects all 3 edges
template<>
template<>
void object::test<7> ()
{
    constexpr int size = 3;
    GEOSGeometry* geoms[size];
    geoms[0] = GEOSGeomFromWKT("LINESTRING(0 0, 10 0)");
    geoms[1] = GEOSGeomFromWKT("LINESTRING(10 0, 20 0)");
    geoms[2] = GEOSGeomFromWKT("LINESTRING(20 0, 30 0)");

    size_t* result = GEOSLongestShortestPath(geoms, size);

    ensure(nullptr != result);
    int count = 0;
    for (int i = 0; i < size; ++i)
        if (result[i] > 0) count++;
    ensure_equals(count, 3);

    GEOSFree(result);
    for (auto& g : geoms) GEOSGeom_destroy(g);
}

// LongestShortestPath: star graph selects 2 edges (one arm through center)
template<>
template<>
void object::test<8> ()
{
    constexpr int size = 4;
    GEOSGeometry* geoms[size];
    geoms[0] = GEOSGeomFromWKT("LINESTRING(0 0, 10 0)");
    geoms[1] = GEOSGeomFromWKT("LINESTRING(0 0, 0 10)");
    geoms[2] = GEOSGeomFromWKT("LINESTRING(0 0, -10 0)");
    geoms[3] = GEOSGeomFromWKT("LINESTRING(0 0, 0 -10)");

    size_t* result = GEOSLongestShortestPath(geoms, size);

    ensure(nullptr != result);
    int count = 0;
    for (int i = 0; i < size; ++i)
        if (result[i] > 0) count++;
    ensure_equals(count, 2);

    GEOSFree(result);
    for (auto& g : geoms) GEOSGeom_destroy(g);
}

// LongestShortestPath: empty input returns NULL
template<>
template<>
void object::test<9> ()
{
    size_t* result = GEOSLongestShortestPath(nullptr, 0);
    ensure(nullptr == result);
}

} // namespace tut
