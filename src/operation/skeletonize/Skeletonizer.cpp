/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2025 Paul Ramsey <pramsey@cleverelephant.ca>
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#include <geos/operation/skeletonize/Skeletonizer.h>
#include <geos/operation/skeletonize/SegmentGraph.h>
#include <geos/operation/skeletonize/CoordinateConditioner.h>
#include <geos/operation/skeletonize/InputOutputs.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/GeometryFilter.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/LineSegment.h>
#include <geos/geom/LineString.h>
#include <geos/geom/MultiLineString.h>
#include <geos/geom/MultiPoint.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/prep/PreparedPolygon.h>
#include <geos/geom/util/Densifier.h>
#include <geos/triangulate/VoronoiDiagramBuilder.h>
#include <geos/util/GEOSException.h>
#include <geos/operation/distance/IndexedFacetDistance.h>
#include <geos/operation/distance/GeometryLocation.h>


#include <cmath>


using namespace geos::geom;
using namespace geos::geom::prep;
using geos::triangulate::VoronoiDiagramBuilder;
using geos::operation::distance::IndexedFacetDistance;
using geos::operation::distance::GeometryLocation;


namespace geos {
namespace operation {   // geos.operation
namespace skeletonize { // geos.operation.skeletonize

Skeletonizer::Skeletonizer(const Polygon &poly, const MultiPoint* points)
        : inputPolygon(poly)
        , inputPoints(points)
        {};

Skeletonizer::Skeletonizer(const Polygon &poly)
        : inputPolygon(poly)
        , inputPoints(nullptr)
        {};


/* public static */
std::unique_ptr<MultiLineString>
Skeletonizer::skeletonize(const Polygon& poly, const MultiPoint& points)
{
    Skeletonizer skel(poly, &points);
    return skel.skeletonize();
}

/* public static */
std::unique_ptr<MultiLineString>
Skeletonizer::skeletonize(const Polygon& poly)
{
    Skeletonizer skel(poly, nullptr);
    return skel.skeletonize();
}

/* private */
bool
Skeletonizer::hasInOutPoints() const
{
        return inputPoints && ! inputPoints->isEmpty();
};


void
Skeletonizer::calculateStatistics(
    const CoordinateSequence* cs,
    SegmentStatistics& stats) const
{
    double n = 0, mean = 0, M2 = 0;
    double max = std::numeric_limits<double>::min();
    double min = std::numeric_limits<double>::max();

    for (std::size_t i = 1; i < cs->size(); i++) {
        const auto& c0 = cs->getAt<CoordinateXY>(i-1);
        const auto& c1 = cs->getAt<CoordinateXY>(i);
        double d = c0.distance(c1);
        max = d > max ? d : max;
        min = d < min ? d : min;
        n += 1;
        double delta = d - mean;
        mean += delta / n;
        M2 += delta * (d - mean);
    }

    stats.numSegments = n;
    stats.averageLength = mean;
    stats.stdevLength = std::sqrt(M2/(n-1));
    stats.minLength = min;
    stats.maxLength = max;
}


void
Skeletonizer::calculateStatistics(
    const Polygon& poly,
    SegmentStatistics& stats) const
{
    const LinearRing* ering = poly.getExteriorRing();
    calculateStatistics(ering->getCoordinatesRO(), stats);

    for (std::size_t i = 0; i < poly.getNumInteriorRing(); i++) {
        const LinearRing* iring = poly.getInteriorRingN(i);
        calculateStatistics(iring->getCoordinatesRO(), stats);
    }
    stats.height = poly.getEnvelopeInternal()->getHeight();
    stats.width = poly.getEnvelopeInternal()->getWidth();
}


std::vector<const Geometry*>
Skeletonizer::findContainedEdges(
    const MultiLineString& allEdges) const
{
    PreparedPolygon preparedInput(&inputPolygon);

    std::vector<const Geometry*> containedEdges;

    struct EdgeFilter : public GeometryFilter {

        PreparedPolygon& filterPoly;
        std::vector<const Geometry*>& edges;

        EdgeFilter(PreparedPolygon& fp, std::vector<const Geometry*>& e)
            : filterPoly(fp)
            , edges(e) {};

        void filter_ro(const Geometry* geom) override {
            if (filterPoly.contains(geom))
                edges.push_back(geom);
        }
    };

    EdgeFilter ef(preparedInput, containedEdges);
    allEdges.apply_ro(&ef);
    return containedEdges;
}


std::vector<GeometryLocation>
Skeletonizer::findInputOutputEdges(
    const MultiLineString& allEdges) const
{
    std::vector<GeometryLocation> inoutEdges;
    if (!hasInOutPoints())
        return inoutEdges;

    //
    // For each input/output point, find the voronoi edge
    // it is closest to. This should be the edge that passes
    // through the gapped points we set up during conditioning.
    //
    std::cout << "findInputOutputEdges " << std::endl;
    std::cout << "  iterating on inputPoints" << std::endl;

    //
    // Index the edges to make the searching faster
    //
    IndexedFacetDistance ifd(&allEdges);

    for (std::size_t i = 0; i < inputPoints->getNumGeometries(); i++) {
        std::cout << "  i = " << i << std::endl;
        const Point* pt = inputPoints->getGeometryN(i);
        std::vector<GeometryLocation> locs = ifd.nearestLocations(pt);

        std::cout << "    locs[0] = " << locs[0].toString() << std::endl;
        std::cout << "    locs[1] = " << locs[1].toString() << std::endl;

        // Because we fed the IndexedFacetDistance a MultLinestring
        // of two-point lines, the geometry component will just be the
        // two-point line we are interested in.
        const Geometry* geomComp = locs[0].getGeometryComponent();
        const std::size_t segIndex = locs[0].getSegmentIndex();
        CoordinateXY inoutCoord = pt->getCoordinatesRO()->getAt(0);
        inoutEdges.emplace_back(geomComp, segIndex, inoutCoord);
    }

    return inoutEdges;
}


/* public */
std::unique_ptr<MultiLineString>
Skeletonizer::skeletonize()
{
    const GeometryFactory* inputFactory = inputPolygon.getFactory();

    //
    // Figure out what we know about this geometry
    //
    SegmentStatistics stats;
    calculateStatistics(inputPolygon, stats);
    double tolerance = 1.0;
    double maxLen = 10.0;

    std::cout << "inputPolygon" << std::endl << inputPolygon << std::endl << std::endl;

    //
    // Condition the input vertices to make the subsequent
    // voronoi results more attractive and less likely to
    // have boundary interactions. Removing very short lines
    // and very acute angles should help.
    //
    std::unique_ptr<MultiLineString> conditionedInput =
        CoordinateConditioner::condition(
            inputPolygon, tolerance);

    std::cout << "conditionedInputs" << std::endl << conditionedInput->toString() << std::endl << std::endl;

    //
    // If user has provided input/output points, we need to
    // replace vertices with gapped points at those locations
    //
    if (hasInOutPoints()) {
        conditionedInput = InputOutputs::addInputOutputGaps(
            *conditionedInput,
            *inputPoints,
            tolerance);

        std::cout << "inputPoints = " << inputPoints->getNumGeometries() << std::endl;
        std::cout << "gappedInputs" << std::endl << conditionedInput->toString() << std::endl << std::endl;
    }

    //
    // Densify the input vertices to make sure that the longest
    // edge length is quite a bit smaller than the overall object
    // size. This makes the voronoi edges, particular the central
    // "skeletal" ones shorter and formed into nice curves more.
    //
    std::unique_ptr<MultiLineString> densifiedInput =
        CoordinateConditioner::densify(
            *conditionedInput, maxLen);

    std::cout << "densifiedInput" << std::endl << densifiedInput->toString() << std::endl << std::endl;

    //
    // This is a point voronoi, and the edges created
    // partition the space so they are equidistant to the
    // input points. This makes some of the edges (the ones
    // near the middle of the shape) excellent candidates to
    // form into a skeleton line.
    //
    VoronoiDiagramBuilder builder;
    //builder.setTolerance(tolerance);
    builder.setSites(*densifiedInput);

    std::unique_ptr<MultiLineString> allVoronoiEdges =
        builder.getDiagramEdges(*inputFactory);

    std::cout << "allVoronoiEdges" << std::endl;
    std::cout << *allVoronoiEdges << std::endl << std::endl;

    //
    // The voronoi edges fully internal to the input polygon
    // form the main body of the skeleton. From these edges
    // we can form a graph and then work out shortest
    // paths through the graph from input to output points
    // and vice versa.
    //
    std::vector<const Geometry*> containedEdges = findContainedEdges(*allVoronoiEdges);

    std::vector<GeometryLocation> inoutEdges = findInputOutputEdges(*allVoronoiEdges);

    auto edgesFiltered = inputFactory->createMultiLineString(containedEdges);
    std::cout << "filteredVoronoiEdges" << std::endl;
    std::cout << *edgesFiltered << std::endl;

    SegmentGraph sg(containedEdges, inoutEdges, inputFactory);
    if (hasInOutPoints()) {
        return sg.shortestPathSkeleton();
    }
    else {
        return sg.longestPathSkeleton();
    }
}



} // namespace geos.operation.skeletonize
} // namespace geos.operation
} // namespace geos
