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
#include <geos/geom/util/Densifier.h>
#include <geos/triangulate/VoronoiDiagramBuilder.h>
#include <geos/util/GEOSException.h>
#include <geos/operation/distance/IndexedFacetDistance.h>
#include <geos/operation/distance/GeometryLocation.h>
#include <geos/simplify/DouglasPeuckerSimplifier.h>

#include <cmath>

#undef DEBUG_OUTPUT

using namespace geos::geom;
using namespace geos::geom::prep;
using geos::triangulate::VoronoiDiagramBuilder;
using geos::operation::distance::IndexedFacetDistance;
using geos::operation::distance::GeometryLocation;
using geos::simplify::DouglasPeuckerSimplifier;


namespace geos {
namespace operation {   // geos.operation
namespace skeletonize { // geos.operation.skeletonize

Skeletonizer::Skeletonizer(const Polygon &poly, const MultiPoint* points)
        : m_inputPolygon(poly)
        , m_inputPoints(points)
        , m_inputFactory(poly.getFactory())
        , m_tolerance(0.0)
        , m_conditioningLength(0.0)
        {}

Skeletonizer::Skeletonizer(const Polygon &poly)
        : m_inputPolygon(poly)
        , m_inputPoints(nullptr)
        , m_inputFactory(poly.getFactory())
        , m_tolerance(0.0)
        , m_conditioningLength(0.0)
        {}


/* public static */
std::unique_ptr<MultiLineString>
Skeletonizer::skeletonize(
    const Polygon& poly,
    double tolerance)
{
    Skeletonizer skel(poly, nullptr);
    auto gf = poly.getFactory();
    return gf->createMultiLineString(skel.skeletonize(tolerance, 0.0));
}


std::unique_ptr<MultiLineString>
Skeletonizer::skeletonize(
    const MultiPolygon& mpoly,
    double tolerance)
{
    std::vector<std::unique_ptr<LineString>> lines;
    for (std::size_t i = 0; i < mpoly.getNumGeometries(); i++) {
        const Polygon* poly = mpoly.getGeometryN(i);
        Skeletonizer skel(*poly, nullptr);
        auto skelLines = skel.skeletonize(tolerance, 0.0);
        for (auto& skln : skelLines) {
            lines.emplace_back(skln.release());
        }
    }
    auto gf = mpoly.getFactory();
    return gf->createMultiLineString(std::move(lines));
}


/* public static */
std::unique_ptr<MultiLineString>
Skeletonizer::skeletonize(
    const Polygon& poly,
    const MultiPoint& mpoints,
    double tolerance)
{
    Skeletonizer skel(poly, &mpoints);
    auto skelLines = skel.skeletonize(tolerance, 0.0);
    auto gf = poly.getFactory();
    return gf->createMultiLineString(std::move(skelLines));
}


/* public static */
std::unique_ptr<MultiLineString>
Skeletonizer::skeletonize(
    const MultiPolygon& mpoly,
    const MultiPoint& mpoints,
    double tolerance)
{
    auto gf = mpoly.getFactory();

    std::vector<std::unique_ptr<LineString>> lines;
    for (std::size_t i = 0; i < mpoly.getNumGeometries(); i++) {

        std::vector<const Geometry*> pts;
        const Polygon* poly = mpoly.getGeometryN(i);

        IndexedFacetDistance ifd(poly);
        for (std::size_t j = 0; j < mpoints.getNumGeometries(); j++) {
            const Point* pt = mpoints.getGeometryN(j);
            if (ifd.isWithinDistance(pt, tolerance)) {
                pts.push_back(pt);
            }
        }

        std::unique_ptr<MultiPoint> inoutPts = gf->createMultiPoint(pts);

        Skeletonizer skel(*poly, inoutPts.get());
        auto skelLines = skel.skeletonize(tolerance, 0.0);

        for (auto& skln : skelLines) {
            lines.emplace_back(skln.release());
        }
    }
    return gf->createMultiLineString(std::move(lines));
}


/* private */
bool
Skeletonizer::hasInOutPoints() const
{
        return m_inputPoints && ! m_inputPoints->isEmpty();
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
    PreparedPolygon preparedInput(&m_inputPolygon);

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
#ifdef DEBUG_OUTPUT
    std::cout << "findInputOutputEdges " << std::endl;
    std::cout << "  iterating on inputPoints" << std::endl;
#endif
    //
    // Index the edges to make the searching faster
    //
    IndexedFacetDistance ifd(&allEdges);

    for (std::size_t i = 0; i < m_inputPoints->getNumGeometries(); i++) {
        const Point* pt = m_inputPoints->getGeometryN(i);
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


/* public */
std::vector<std::unique_ptr<LineString>>
Skeletonizer::skeletonize(double tolerance, double conditioningLength)
{
    setTolerance(tolerance);
    setConditioningLength(conditioningLength);
    return skeletonize();
}


/* public */
std::vector<std::unique_ptr<LineString>>
Skeletonizer::skeletonize()
{
    //
    // Figure out what we know about this geometry
    //
    if (m_tolerance <= 0.0 || m_conditioningLength <= 0) {
        SegmentStatistics stats = calculateStatistics(m_inputPolygon);
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

    //
    // Condition the input vertices to make the subsequent
    // voronoi results more attractive and less likely to
    // have boundary interactions. Removing very short lines
    // and very acute angles should help.
    //
    std::unique_ptr<MultiLineString> conditionedInput =
        CoordinateConditioner::condition(
            m_inputPolygon, m_tolerance);

#ifdef DEBUG_OUTPUT
    std::cout << "conditionedInput" << std::endl << conditionedInput->toString() << std::endl << std::endl;
#endif

    //
    // If user has provided input/output points, we need to
    // replace vertices with gapped points at those locations
    //
    if (hasInOutPoints()) {
        conditionedInput = InputOutputs::addInputOutputGaps(
            *conditionedInput,
            *m_inputPoints,
            m_tolerance);

#ifdef DEBUG_OUTPUT
        std::cout << "m_inputPoints = " << m_inputPoints->getNumGeometries() << std::endl;
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
    std::vector<const Geometry*> containedEdges = findContainedEdges(*allVoronoiEdges);

    std::vector<GeometryLocation> inoutEdges = findInputOutputEdges(*allVoronoiEdges);

#ifdef DEBUG_OUTPUT
    std::cout << "filteredVoronoiEdges" << std::endl;
    std::cout << *(getGeometry(containedEdges)) << std::endl;
    std::cout << "inoutVoronoiEdges" << std::endl;
    std::cout << *(getGeometry(inoutEdges)) << std::endl;
#endif

    SegmentGraph sg(containedEdges, inoutEdges, m_inputFactory);
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
