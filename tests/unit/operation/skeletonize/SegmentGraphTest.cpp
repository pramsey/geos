//
// Test Suite for geos::operation::skeletonize::SegmentGraph class.
//

// tut
#include <tut/tut.hpp>
// geos
#include <geos/operation/skeletonize/SegmentGraph.h>
#include <geos/operation/distance/GeometryLocation.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/MultiLineString.h>
#include <geos/geom/MultiPoint.h>
#include <geos/geom/MultiPolygon.h>
#include <geos/geom/Polygon.h>
#include <geos/io/WKBReader.h>
#include <geos/io/WKTReader.h>
#include <geos/io/WKTWriter.h>
// std
#include <memory>
#include <string>
#include <vector>
#include <iostream>

namespace tut {
//
// Test Group
//

using namespace geos::geom;
using namespace geos::operation::skeletonize;
using geos::operation::distance::GeometryLocation;

// Common data used by tests
struct test_segmentgraph_data {

    geos::io::WKTReader _r;
    geos::io::WKTWriter _w;
    std::unique_ptr<MultiLineString> edgeVectorGeometry;

    test_segmentgraph_data() {};

    std::unique_ptr<MultiPoint>
    read_mpoint(const std::string& wkt)
    {
        std::unique_ptr<MultiPoint> g = _r.read<MultiPoint>(wkt);
        return g;
    }

    std::unique_ptr<Polygon>
    read_polygon(const std::string& wkt)
    {
        std::unique_ptr<Polygon> g = _r.read<Polygon>(wkt);
        return g;
    }

    std::unique_ptr<MultiPolygon>
    read_mpolygon(const std::string& wkt)
    {
        std::unique_ptr<MultiPolygon> g = _r.read<MultiPolygon>(wkt);
        return g;
    }

    std::unique_ptr<MultiLineString>
    read_mline(const std::string& wkt)
    {
        std::unique_ptr<MultiLineString> g = _r.read<MultiLineString>(wkt);
        return g;
    }

    std::vector<const Geometry*>
    read_edge_vector(const std::string& wkt)
    {
        edgeVectorGeometry = read_mline(wkt);
        std::vector<const Geometry*> edges;
        for (std::size_t i = 0; i < edgeVectorGeometry->getNumGeometries(); i++) {
            edges.push_back(edgeVectorGeometry->getGeometryN(i));
        }
        return edges;
    }

};

typedef test_group<test_segmentgraph_data> group;
typedef group::object object;

group test_segmentgraph_group("geos::operation::skeletonize::SegmentGraph");


template<>
template<>
void object::test<1> ()
{
    auto edges = read_edge_vector(
        "MULTILINESTRING ((0 0, 1 1), (1 1, 2 2), (1 1, 2 1))");

    std::vector<GeometryLocation> locs;
    auto gf = geos::geom::GeometryFactory::create();

    SegmentGraph sg(edges, locs, gf.get());
    sg.build();

    auto endVertices = sg.endVertices();
    ensure_equals("endVertices == 3", endVertices.size(), 3ul);

    auto result1 = sg.shortestPath(endVertices[0], endVertices[1]);
    ensure_equals("shortestPath.second", result1.second, 2.41421, 0.0001);

    auto result2 = sg.longestPath(endVertices);
    ensure_equals("longestPath.second", result2.second, 2.82843, 0.0001);
}

template<>
template<>
void object::test<2>()
{
    // Two separate line segments, disconnected
    auto edges = read_edge_vector(
        "MULTILINESTRING ((0 0, 1 0), (10 0, 11 0))");

    std::vector<GeometryLocation> locs;
    auto gf = geos::geom::GeometryFactory::create();

    SegmentGraph sg(edges, locs, gf.get());
    sg.build();

    const auto& sgConnected = sg.connectedComponents();
    uint32_t sgConnectedCount = sg.connectedComponentsCount();

    ensure_equals("number of connected components", sgConnected.size(), 2ul);

    // Each connected component should have exactly 2 vertices
    for (const auto& comp : sgConnected) {
        ensure_equals("component vertex count", comp.second.size(), 2ul);
    }

    // m_connectedCount should match m_connected.size()
    ensure_equals("m_connectedCount", sgConnectedCount, sgConnected.size());
}

template<>
template<>
void object::test<3>()
{
    // Empty graph: no segments
    std::vector<const Geometry*> edges;
    std::vector<GeometryLocation> locs;
    auto gf = geos::geom::GeometryFactory::create();

    SegmentGraph sg(edges, locs, gf.get());
    sg.build();

    ensure_equals("no components in empty graph", sg.connectedComponents().size(), 0ul);
    ensure_equals("m_connectedCount == 0", sg.connectedComponentsCount(), 0u);
}

template<>
template<>
void object::test<4>()
{
    // Single edge: should result in one component of two vertices
    auto edges = read_edge_vector("MULTILINESTRING ((0 0, 1 0))");
    std::vector<GeometryLocation> locs;
    auto gf = geos::geom::GeometryFactory::create();

    SegmentGraph sg(edges, locs, gf.get());
    sg.build();

    const auto& comps = sg.connectedComponents();
    ensure_equals("one component", comps.size(), 1ul);
    ensure_equals("m_connectedCount == 1", sg.connectedComponentsCount(), 1u);

    for (const auto& c : comps) {
        ensure_equals("two vertices in component", c.second.size(), 2ul);
    }
}

template<>
template<>
void object::test<5>()
{
    // Connected loop: triangle, all vertices part of one component
    auto edges = read_edge_vector(
        "MULTILINESTRING ((0 0, 1 0), (1 0, 0.5 1), (0.5 1, 0 0))"
    );
    std::vector<GeometryLocation> locs;
    auto gf = geos::geom::GeometryFactory::create();

    SegmentGraph sg(edges, locs, gf.get());
    sg.build();

    const auto& comps = sg.connectedComponents();
    uint32_t sgConnectedCount = sg.connectedComponentsCount();
    ensure_equals("one triangle component", comps.size(), 1ul);
    ensure_equals("m_connectedCount == 1", sgConnectedCount, 1u);

    for (const auto& c : comps) {
        ensure_equals("three vertices in triangle", c.second.size(), 3ul);
    }
}



} // namespace tut

