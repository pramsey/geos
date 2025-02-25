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

#include <cmath>


using namespace geos::geom;
using namespace geos::geom::prep;
using geos::triangulate::VoronoiDiagramBuilder;


namespace geos {
namespace operation {   // geos.operation
namespace skeletonize { // geos.operation.skeletonize

Skeletonizer::Skeletonizer(const Polygon &poly, const MultiPoint& points)
        : inputPolygon(poly)
        , inputPoints(points)
        {};


/* public static */
std::unique_ptr<LineString>
Skeletonizer::skeletonize(const Polygon& poly, const MultiPoint& points)
{
    Skeletonizer skel(poly, points);
    return skel.skeletonize();
}


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


/* public */
std::unique_ptr<LineString>
Skeletonizer::skeletonize()
{
    const GeometryFactory* inputFactory = inputPolygon.getFactory();

    // Figure out what we know about this geometry
    SegmentStatistics stats;
    calculateStatistics(inputPolygon, stats);
    double tolerance = 1.0;
    double maxLen = 10.0;

    std::cout << "inputPolygon" << std::endl << inputPolygon << std::endl << std::endl;

    // Condition the input vertices by
    // removing spikes/gores
    // removing dupes and very short segments
    // adding extra segments to force maximum segment length
    std::unique_ptr<Geometry> conditionedInput =
        CoordinateConditioner::condition(
            &inputPolygon, tolerance, maxLen);

    std::cout << "conditionedInput" << std::endl << *conditionedInput << std::endl << std::endl;

    VoronoiDiagramBuilder builder;
    //builder.setTolerance(tolerance);
    builder.setSites(*conditionedInput);

    std::unique_ptr<MultiLineString> allEdges = builder.getDiagramEdges(*inputFactory);

    std::cout << "allEdges" << std::endl;
    std::cout << *allEdges << std::endl << std::endl;

    std::vector<const Geometry*> nonCrossingEdges;
    PreparedPolygon preparedInput(&inputPolygon);

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

    EdgeFilter ef(preparedInput, nonCrossingEdges);
    allEdges->apply_ro(&ef);

    auto edgesFiltered = inputFactory->createMultiLineString(nonCrossingEdges);
    std::cout << "edgesFiltered" << std::endl;
    std::cout << *edgesFiltered << std::endl;

    if (inputPoints.getNumGeometries() > 0)
        std::cout << "inputPoints.getNumGeometries() > 0" << std::endl;

    SegmentGraph sg(*edgesFiltered);
    auto skeletonLine = sg.longestPath();

    return skeletonLine;
}



} // namespace geos.operation.skeletonize
} // namespace geos.operation
} // namespace geos
