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
#include <geos/algorithm/Angle.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/GeometryFilter.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/LineSegment.h>
#include <geos/geom/LineString.h>
#include <geos/geom/MultiLineString.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/prep/PreparedPolygon.h>
#include <geos/geom/util/Densifier.h>
#include <geos/triangulate/VoronoiDiagramBuilder.h>
#include <geos/util/GEOSException.h>

using geos::algorithm::Angle;
using geos::geom::Coordinate;
using geos::geom::CoordinateXY;
using geos::geom::CoordinateSequence;
using geos::geom::GeometryTypeId;
using geos::geom::Geometry;
using geos::geom::GeometryFactory;
using geos::geom::GeometryFilter;
using geos::geom::LinearRing;
using geos::geom::LineSegment;
using geos::geom::LineString;
using geos::geom::MultiLineString;
using geos::geom::Polygon;
using geos::geom::prep::PreparedPolygon;
using geos::geom::util::Densifier;
using geos::triangulate::VoronoiDiagramBuilder;


namespace geos {
namespace operation {   // geos.operation
namespace skeletonize { // geos.operation.skeletonize


/* public static */
std::unique_ptr<LineString>
Skeletonizer::skeletonize(const Polygon& poly)
{
    Skeletonizer skel(poly);
    return skel.skeletonize();
}


double
Skeletonizer::startPointAngle(const CoordinateSequence* cs) const
{
    CoordinateXY p0, p1, p2;
    cs->getAt(cs->size()-2, p0); // start/end point are dupes
    cs->getAt(0, p1);
    cs->getAt(1, p2);
    return Angle::angleBetween(p0, p1, p2);
}


std::size_t
Skeletonizer::findWidestAngle(const CoordinateSequence* cs) const
{
    std::size_t widest = 0;
    double widestAngle = 0.0;
    for (std::size_t i = 0; i < cs->size()-1; i++) {
        CoordinateXY p0, p1, p2;
        if (i > 0) {
            cs->getAt(i-1, p0);
            cs->getAt(i,   p1);
            cs->getAt(i+1, p2);
        }
        else {
            cs->getAt(cs->size()-2, p0);
            cs->getAt(0, p1);
            cs->getAt(1, p2);
        }
        double angle = Angle::angleBetween(p0, p1, p2);
        if (angle > widestAngle) {
            widestAngle = angle;
            widest = i;
        }
    }
    return widest;
}


std::unique_ptr<CoordinateSequence>
Skeletonizer::startSequenceAtWidest(const CoordinateSequence* seq)
{
    std::size_t widest = findWidestAngle(seq);
    auto reorientedSeq = std::make_unique<CoordinateSequence>(seq->size());
    // Start at widest and take everything up to the end
    for (std::size_t i = widest; i < seq->size(); i++) {
        reorientedSeq->add(seq->getAt(i));
    }
    // Skip first point (dupes the last) and also
    // take widest point again, so new ring has dupe
    // start/end points
    for (std::size_t i = 1; i <= widest; i++) {
        reorientedSeq->add(seq->getAt(i));
    }
    return reorientedSeq;
}


std::unique_ptr<Geometry>
Skeletonizer::densifyDefault(const Polygon* poly, double tolerance)
{
    return Densifier::densify(poly, tolerance);
}


std::unique_ptr<Geometry>
Skeletonizer::densifyUniformly(const Polygon* poly, double tolerance)
{
    auto denseExtRing = densifyUniformly(poly->getExteriorRing(), tolerance);
    const GeometryFactory* inputFactory = inputPolygon.getFactory();

    std::vector<std::unique_ptr<LinearRing>> denseIntRings;
    for (std::size_t i = 0; i < poly->getNumInteriorRing(); i++) {
        auto denseRing = densifyUniformly(poly->getInteriorRingN(i), tolerance);
        denseIntRings.emplace_back(denseRing.release());
    }
    std::unique_ptr<Polygon> result = inputFactory->createPolygon(std::move(denseExtRing), std::move(denseIntRings));
    return result;
}


std::unique_ptr<LinearRing>
Skeletonizer::densifyUniformly(const LinearRing* ring, double tolerance)
{
    const GeometryFactory* inputFactory = inputPolygon.getFactory();
    const CoordinateSequence* coords = ring->getCoordinatesRO();
    auto denseCoords = std::make_unique<CoordinateSequence>();

    LineSegment seg;
    double remainder = 0.0;
    for (std::size_t i = 1; i < coords->size(); i++) {
        seg.p0 = coords->getAt<Coordinate>(i-1);
        seg.p1 = coords->getAt<Coordinate>(i);
        double segLen = seg.getLength();

        // First coordinate is always added
        denseCoords->add(seg.p0, false);

        // No room to densify this small segment
        if (tolerance >= segLen)
            continue;

        // Add first sub-point
        Coordinate c;
        seg.pointAlong(remainder / segLen, c);
        denseCoords->add(c, false);
        remainder = segLen - remainder;

        // Add as many other points as we can
        while (remainder > tolerance) {
            remainder -= tolerance;
            seg.pointAlong(1.0 - (remainder/segLen), c);
            denseCoords->add(c, false);
        }
    }

    // add final point in ring
    denseCoords->add(seg.p1, false);

    return inputFactory->createLinearRing(std::move(denseCoords));
}

/* public */
std::unique_ptr<LineString>
Skeletonizer::skeletonize()
{
    std::cout << "GeometryTypeId == " << inputPolygon.getGeometryTypeId() << std::endl;
    const GeometryFactory* inputFactory = inputPolygon.getFactory();

    // Naive densifier
    // auto denseInput = densifyDefault(&inputPolygon, 5.0);
    auto denseInput = densifyUniformly(&inputPolygon, 6.0);

    std::cout << *denseInput << std::endl << std::endl;

    VoronoiDiagramBuilder builder;
    //builder.setTolerance(tolerance);
    builder.setSites(*denseInput);

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

    SegmentGraph sg(*edgesFiltered);
    auto skeletonLine = sg.longestPath();

    return skeletonLine;
}



} // namespace geos.operation.skeletonize
} // namespace geos.operation
} // namespace geos
