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
#include <geos/geom/MultiPolygon.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/prep/PreparedPolygon.h>
#include <geos/geom/util/PointExtracter.h>
#include <geos/geom/util/PolygonExtracter.h>
#include <geos/triangulate/VoronoiDiagramBuilder.h>
#include <geos/util/GEOSException.h>
#include <geos/operation/distance/IndexedFacetDistance.h>
#include <geos/operation/distance/GeometryLocation.h>
#include <geos/simplify/DouglasPeuckerSimplifier.h>

#include <cmath>

#undef DEBUG_OUTPUT

using namespace geos::geom;
using namespace geos::geom::prep;
using geos::geom::util::PointExtracter;
using geos::geom::util::PolygonExtracter;
using geos::triangulate::VoronoiDiagramBuilder;
using geos::operation::distance::IndexedFacetDistance;
using geos::operation::distance::GeometryLocation;
using geos::simplify::DouglasPeuckerSimplifier;


namespace geos {
namespace operation {   // geos.operation
namespace skeletonize { // geos.operation.skeletonize

Skeletonizer::Skeletonizer(const Geometry* poly, const Geometry* points)
        : m_inputGeometry(poly)
        , m_inputPoints(points)
        , m_inputFactory(poly->getFactory())
        , m_tolerance(0.0)
        , m_conditioningLength(0.0)
        {}

Skeletonizer::Skeletonizer(const Geometry *poly)
        : m_inputGeometry(poly)
        , m_inputPoints(nullptr)
        , m_inputFactory(poly->getFactory())
        , m_tolerance(0.0)
        , m_conditioningLength(0.0)
        {}


/* public static */
std::unique_ptr<MultiLineString>
Skeletonizer::skeletonize(
    const Geometry* poly,
    const Geometry* mpoints,
    double tolerance,
    double conditioningLength)
{
    Skeletonizer skel(poly, mpoints);
    skel.setTolerance(tolerance);
    skel.setConditioningLength(conditioningLength);
    return skel.skeletonize();
}


/* public static */
std::unique_ptr<MultiLineString>
Skeletonizer::skeletonize(
    const Geometry* poly,
    double tolerance,
    double conditioningLength)
{
    Skeletonizer skel(poly, nullptr);
    skel.setTolerance(tolerance);
    skel.setConditioningLength(conditioningLength);
    return skel.skeletonize();
}


void
Skeletonizer::initializeStatistics(
    SegmentStatistics& stats) const
{
    stats.numSegments = 0;
    stats.maxLength = std::numeric_limits<double>::max();
    stats.minLength = std::numeric_limits<double>::min();
    stats.averageLength = 0.0;
    stats.M2 = 0.0;
    stats.stdevLength = 0.0;
    stats.width = 0.0;
    stats.height = 0.0;
    stats.maxX = std::numeric_limits<double>::min();
    stats.maxY = std::numeric_limits<double>::min();
}


void
Skeletonizer::calculateStatistics(
    const CoordinateSequence* cs,
    SegmentStatistics& stats) const
{
    for (std::size_t i = 1; i < cs->size(); i++) {
        const auto& c0 = cs->getAt<CoordinateXY>(i-1);
        const auto& c1 = cs->getAt<CoordinateXY>(i);
        double len = c0.distance(c1);

        if (len > stats.maxLength) stats.maxLength = len;
        if (len < stats.minLength) stats.minLength = len;

        if (c0.x > stats.maxX) stats.maxX = c0.x;
        if (c0.y > stats.maxY) stats.maxY = c0.y;

        stats.numSegments += 1;
        double delta = len - stats.averageLength;
        stats.averageLength += delta / stats.numSegments;
        stats.M2 += delta * (len - stats.averageLength);
    }
}


Skeletonizer::SegmentStatistics
Skeletonizer::calculateStatistics(
    const Polygon& poly) const
{
    SegmentStatistics stats;
    initializeStatistics(stats);

    const LinearRing* ering = poly.getExteriorRing();
    calculateStatistics(ering->getCoordinatesRO(), stats);

    for (std::size_t i = 0; i < poly.getNumInteriorRing(); i++) {
        const LinearRing* iring = poly.getInteriorRingN(i);
        calculateStatistics(iring->getCoordinatesRO(), stats);
    }
    stats.height = poly.getEnvelopeInternal()->getHeight();
    stats.width = poly.getEnvelopeInternal()->getWidth();
    stats.stdevLength = std::sqrt(stats.M2/(stats.averageLength-1));
    return stats;
}

double
Skeletonizer::defaultTolerance(
    const SegmentStatistics& stats) const
{
    double coordinateMagnitude = std::fabs((stats.maxX + stats.maxY) / 2.0);
    double absMag = std::fabs(coordinateMagnitude);
    // Probably geographics, so tolerance of about 10cm
    if (absMag < 360)
        return 0.0000001;
    // Planar so tolerance of 1/10 of a unit (meter, foot) seems
    // reasonable
    else
        return 0.1;
}

double
Skeletonizer::defaultConditioningLength(
    const SegmentStatistics& stats) const
{
    return std::fabs(stats.averageLength - stats.stdevLength) / 10;
}


//
// Given all the edges generated by the Voronoi,
// extract out just those edges that are fully
// contained by the input polygon. These will
// be used to form the skeleton.
//
std::vector<const Geometry*>
Skeletonizer::findContainedEdges(
    const MultiLineString& allEdges) const
{
    PreparedPolygon preparedInput(m_inputGeometry);

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

//
// We start with input/output points, but we need to
// associate those points with edges from the Voronoi.
// These edges will guide the joining of the internal
// skeleton edges to the input/output points. The
// input/output edges are, for each input/output point,
// the edge closest to that point.
//
std::vector<GeometryLocation>
Skeletonizer::findInputOutputEdges(
    const MultiLineString& allEdges,
    const std::vector<const Point*>& points) const
{
    std::vector<GeometryLocation> inoutEdges;
    if (points.empty())
        return inoutEdges;

    //
    // For each input/output point, find the voronoi edge
    // it is closest to. This should be the edge that passes
    // through the gapped points we set up during conditioning.
    //
#ifdef DEBUG_OUTPUT
    std::cout << "findInputOutputEdges " << std::endl;
    std::cout << "  iterating on inputPoints" << std::endl;
#endif
    //
    // Index the edges to make the searching faster
    //
    IndexedFacetDistance ifd(&allEdges);

    for (const Point* pt : points) {
        std::vector<GeometryLocation> locs = ifd.nearestLocations(pt);

#ifdef DEBUG_OUTPUT
        std::cout << "  i = " << i << std::endl;
        std::cout << "    locs[0] = " << locs[0].toString() << std::endl;
        std::cout << "    locs[1] = " << locs[1].toString() << std::endl;
#endif
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

//
// Convert a vector of GeometryLocation to an equivalent
// MultiLineString
//
std::unique_ptr<MultiLineString>
Skeletonizer::getGeometry(
    std::vector<GeometryLocation>& inoutEdges) const
{
    std::vector<std::unique_ptr<Geometry>> edges;
    for (const GeometryLocation& gl : inoutEdges) {
        edges.emplace_back(gl.getGeometryComponent()->clone().release());
    }
    return m_inputFactory->createMultiLineString(std::move(edges));
}

//
// Convert a vector of pointers to edges to an equivalent
// MultiLineString
//
std::unique_ptr<MultiLineString>
Skeletonizer::getGeometry(
    std::vector<const Geometry*>& edges) const
{
    return m_inputFactory->createMultiLineString(edges);
}




std::unique_ptr<MultiLineString>
Skeletonizer::skeletonize()
{
    std::vector<const Polygon*> polys;
    PolygonExtracter pe(polys);
    m_inputGeometry->apply_ro(&pe);

    std::vector<const Point*> inputPoints;
    if (m_inputPoints && ! m_inputPoints->isEmpty()) {
        PointExtracter pte(inputPoints);
        m_inputPoints->apply_ro(&pte);
    }

    if (polys.empty()) {
        throw geos::util::GEOSException("No polygons in input geometry");
    }

    std::vector<std::unique_ptr<LineString>> lines;
    for (const Polygon* poly : polys) {

        std::vector<const Point*> points;
        if (!inputPoints.empty()) {
            IndexedFacetDistance ifd(poly);
            for (const Point* pt : inputPoints) {
                if (ifd.isWithinDistance(pt, getTolerance())) {
                    points.push_back(pt);
                }
            }
        }
        skeletonizePolygon(poly, points, lines);
    }
    return m_inputFactory->createMultiLineString(std::move(lines));
}


/* private */
void
Skeletonizer::skeletonizePolygon(
    const Polygon* poly,
    const std::vector<const Point*>& points,
    std::vector<std::unique_ptr<LineString>>& skelLines)
{
    //
    // Figure out what we know about this geometry
    //
    if (m_tolerance <= 0.0 || m_conditioningLength <= 0) {
        SegmentStatistics stats = calculateStatistics(*poly);
        if (m_tolerance <= 0.0)
            m_tolerance = defaultTolerance(stats);

        if (m_conditioningLength <= 0.0)
            m_conditioningLength = defaultConditioningLength(stats);
    }

#ifdef DEBUG_OUTPUT
    std::cout << "m_tolerance = " << m_tolerance << std::endl;
    std::cout << "m_conditioningLength = " << m_conditioningLength << std::endl;
    std::cout << "inputPolygon" << std::endl << m_inputPolygon << std::endl << std::endl;
#endif

    // std::cout << "m_inputGeometry->getNumPoints() = " << m_inputGeometry->getNumPoints() << std::endl;

    //
    // Condition the input vertices to make the subsequent
    // voronoi results more attractive and less likely to
    // have boundary interactions. Removing very short lines
    // and very acute angles should help.
    //
    std::unique_ptr<MultiLineString> conditionedInput =
        CoordinateConditioner::condition(
            *poly, m_tolerance);

#ifdef DEBUG_OUTPUT
    std::cout << "conditionedInput" << std::endl << conditionedInput->toString() << std::endl << std::endl;
#endif

    // std::cout << "conditionedInput->getNumPoints() = " << conditionedInput->getNumPoints() << std::endl;

    //
    // If user has provided input/output points, we need to
    // replace vertices with gapped points at those locations
    //
    if (! points.empty()) {
        conditionedInput = InputOutputs::addInputOutputGaps(
            *conditionedInput,
            points,
            m_tolerance);

#ifdef DEBUG_OUTPUT
        std::cout << "points = " << points->size() << std::endl;
        std::cout << "gappedConditionedInput" << std::endl << conditionedInput->toString() << std::endl << std::endl;
#endif
    }

    //
    // Densify the input vertices to make sure that the longest
    // edge length is quite a bit smaller than the overall object
    // size. This makes the voronoi edges, particular the central
    // "skeletal" ones shorter and formed into nice curves more.
    //
    std::unique_ptr<MultiLineString> densifiedInput =
        CoordinateConditioner::densify(
            *conditionedInput, m_conditioningLength);

#ifdef DEBUG_OUTPUT
    std::cout << "densifiedInput" << std::endl << densifiedInput->toString() << std::endl << std::endl;
#endif

    // std::cout << "densifiedInput->getNumPoints() = " << densifiedInput->getNumPoints() << std::endl;

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
        builder.getDiagramEdges(*m_inputFactory);

#ifdef DEBUG_OUTPUT
    std::cout << "allVoronoiEdges" << std::endl;
    std::cout << *allVoronoiEdges << std::endl << std::endl;
#endif

    //
    // The voronoi edges fully internal to the input polygon
    // form the main body of the skeleton. From these edges
    // we can form a graph and then work out shortest
    // paths through the graph from input to output points
    // and vice versa.
    //
    std::vector<const Geometry*> containedEdges = findContainedEdges(
        *allVoronoiEdges);

    std::vector<GeometryLocation> inoutEdges = findInputOutputEdges(
        *allVoronoiEdges,
        points);

#ifdef DEBUG_OUTPUT
    std::cout << "filteredVoronoiEdges" << std::endl;
    std::cout << *(getGeometry(containedEdges)) << std::endl;
    std::cout << "inoutVoronoiEdges" << std::endl;
    std::cout << *(getGeometry(inoutEdges)) << std::endl;
#endif

    //
    // If we have no input/output points, just return
    // longest path, otherwise find the shortest paths
    // that connect all the input/output points.
    //
    SegmentGraph sg(containedEdges, inoutEdges, m_inputFactory);
    std::vector<std::unique_ptr<LineString>> polySkelLines
        = points.empty()
        ? sg.longestPathSkeleton()
        : sg.shortestPathSkeleton();

    for (auto& ln : polySkelLines) {
        skelLines.emplace_back(ln.release());
    }

    return;
}



} // namespace geos.operation.skeletonize
} // namespace geos.operation
} // namespace geos
