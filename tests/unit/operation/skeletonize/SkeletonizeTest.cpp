//
// Test Suite for geos::operation::skeletonize::Skeletonizer class.
//

// tut
#include <tut/tut.hpp>
// geos
#include <geos/operation/skeletonize/Skeletonizer.h>
#include <geos/operation/skeletonize/SegmentGraph.h>
#include <geos/operation/skeletonize/InputOutputs.h>
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
struct test_skeletonizetest_data {

    geos::io::WKTReader _r;
    geos::io::WKTWriter _w;
    std::unique_ptr<MultiLineString> edgeVectorGeometry;

    test_skeletonizetest_data() {};

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

typedef test_group<test_skeletonizetest_data> group;
typedef group::object object;

group test_skeletonizetest_group("geos::operation::skeletonize::skeletonize");



template<>
template<>
void object::test<1> ()
{
    auto poly = read_polygon("POLYGON ((100 100, 200 100, 200 200, 100 200, 100 100))");
    auto mpoint = read_mpoint("MULTIPOINT EMPTY");
    auto output = Skeletonizer::skeletonize(*poly, *mpoint);
    ensure("Skeleton Positive Length", output->getLength() > 100);
}


template<>
template<>
void object::test<2> ()
{
    auto poly = read_polygon("POLYGON ((100 100, 200 100, 200 200, 100 200, 100 100))");
    auto mpoint = read_mpoint("MULTIPOINT ((200 200))");
    auto output = Skeletonizer::skeletonize(*poly, *mpoint);
    // std::cout << _w.write(*output) << std::endl;
    ensure("Skeleton Positive Length", output->getLength() > 100);
}

template<>
template<>
void object::test<3> ()
{
    auto poly = read_polygon("POLYGON ((100 100, 200 100, 200 200, 100 200, 100 100))");
    auto mpoint = read_mpoint("MULTIPOINT ((200 200), (200 100))");
    auto output = Skeletonizer::skeletonize(*poly, *mpoint);
    // std::cout << _w.write(*output) << std::endl;
    ensure("Skeleton Positive Length", output->getLength() > 100);
}

template<>
template<>
void object::test<4> ()
{
    auto poly = read_polygon("POLYGON ((100 100, 200 100, 200 200, 100 200, 100 100))");
    auto mpoint = read_mpoint("MULTIPOINT (100 100, 200 100, 200 200, 100 200)");
    auto output = Skeletonizer::skeletonize(*poly, *mpoint);
    // std::cout << _w.write(*output) << std::endl;
    ensure("Skeleton Positive Length", output->getLength() > 200);
}

template<>
template<>
void object::test<5> ()
{
    auto poly = read_polygon("POLYGON ((100 100, 200 100, 200 200, 150 200, 100 200, 100 100))");
    auto mpoint = read_mpoint("MULTIPOINT (100 100, 200 100, 155 200)");
    auto output = Skeletonizer::skeletonize(*poly, *mpoint);
    // std::cout << _w.write(*output) << std::endl;
    ensure("Skeleton Positive Length", output->getLength() > 100);
}


template<>
template<>
void object::test<11> ()
{
    auto poly = read_polygon("POLYGON ((100 130, 120 160, 160 180, 270 180, 330 180, 420 170, 430 120, 450 80, 430 50, 350 50, 260 30, 170 30, 130 30, 90 80, 100 130))");
    auto mpoint = read_mpoint("MULTIPOINT((450 80),(160 180), (100 130))");
    auto output = Skeletonizer::skeletonize(*poly, *mpoint);
    // std::cout << _w.write(*output) << std::endl;
}


template<>
template<>
void object::test<12> ()
{
    auto edges = read_edge_vector("MULTILINESTRING ((0 0, 1 1), (1 1, 2 2), (1 1, 2 1))");
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
void object::test<13> ()
{
    auto poly = read_mpolygon("MULTIPOLYGON (((100 130, 120 160, 160 180, 270 180, 330 180, 420 170, 430 120, 450 80, 430 50, 350 50, 260 30, 170 30, 130 30, 90 80, 100 130)))");
    auto mpoint = read_mpoint("MULTIPOINT((450 80),(160 180), (100 130))");
    auto output = Skeletonizer::skeletonize(*poly, *mpoint);
    // std::cout << _w.write(*output) << std::endl;
}


template<>
template<>
void object::test<15> ()
{
    auto poly = read_polygon("POLYGON ((0 0, 100 0, 100 20, 0 20, 0 0))");
    auto pts = read_mpoint("MULTIPOINT (0 0, 100.0 20, 100 0, 0 20)");
    auto output = Skeletonizer::skeletonize(*poly, *pts);
    // std::cout << std::endl << "Skeleton" << std::endl << output->toString() << std::endl;
}



} // namespace tut

