//
// Test Suite for geos::operation::skeletonize::Skeletonizer class.
//

// tut
#include <tut/tut.hpp>
// geos
#include <geos/operation/skeletonize/Skeletonizer.h>
#include <geos/operation/skeletonize/SegmentGraph.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/MultiLineString.h>
#include <geos/geom/Polygon.h>
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

using geos::geom::Geometry;
using geos::geom::Polygon;
using geos::geom::MultiLineString;

// Common data used by tests
struct test_skeletonizetest_data {

    geos::io::WKTReader wktreader;
    geos::io::WKTWriter wktwriter;

    test_skeletonizetest_data()
        : wktreader()
    {
    }

};

typedef test_group<test_skeletonizetest_data> group;
typedef group::object object;

group test_skeletonizetest_group("geos::operation::skeletonize::skeletonize");

template<>
template<>
void object::test<1> ()
{
    // auto geom = wktreader.read("POLYGON ((100 130, 120 160, 160 180, 270 180, 330 180, 420 170, 430 120, 450 80, 430 50, 350 50, 260 30, 170 30, 130 30, 90 80, 100 130))");
    auto geom = wktreader.read("POLYGON ((100 130, 120 160, 160 180, 270 180, 330 180, 420 170, 430 120, 450 80, 430 50, 350 50, 260 30, 170 30, 130 30, 90 80, 100 130), (200 130, 170 80, 250 80, 370 90, 320 140, 200 130))");
    // auto geom = wktreader.read("POLYGON ((120 160, 120.5 160.5, 160 180, 100.05 129.95, 100 130, 120 160))");
    // auto geom = wktreader.read("POLYGON ((100 100, 200 100, 200 200, 100 200, 100 100), (140 140, 160 140, 160 160, 140 160, 140 140))");
    auto poly = static_cast<const Polygon*>(geom.get());
    geos::operation::skeletonize::Skeletonizer skel(*poly);
    auto output = skel.skeletonize();
    std::cout << *output << std::endl;
}


template<>
template<>
void object::test<2> ()
{
    auto geom = wktreader.read("MULTILINESTRING ((0 0, 1 1), (1 1, 2 2), (1 1, 2 1))");
    auto mline = static_cast<const MultiLineString*>(geom.get());

    geos::operation::skeletonize::SegmentGraph sg(*mline);
    sg.buildAdjacencyList();
    std::cout << sg << std::endl;
}


} // namespace tut

