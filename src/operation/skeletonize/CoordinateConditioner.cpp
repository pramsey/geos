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

#include <geos/operation/skeletonize/CoordinateConditioner.h>
#include <geos/algorithm/Angle.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/LineString.h>
#include <geos/geom/MultiLineString.h>
#include <geos/geom/Polygon.h>

#include <queue>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using geos::algorithm::Angle;
using geos::geom::CoordinateXY;
using geos::geom::CoordinateSequence;
using geos::geom::Geometry;
using geos::geom::GeometryFactory;
using geos::geom::Polygon;
using geos::geom::LinearRing;
using geos::geom::LineString;
using geos::geom::MultiLineString;


namespace geos {
namespace operation {   // geos.operation
namespace skeletonize { // geos.operation.skeletonize


/* public static */
std::unique_ptr<MultiLineString>
CoordinateConditioner::condition(
    const Polygon& poly,
    double tolerance)
{
    CoordinateConditioner cc;
    return cc.conditionPolygon(poly, tolerance);
}


/* private */
std::unique_ptr<MultiLineString>
CoordinateConditioner::conditionPolygon(
    const Polygon& poly,
    double tolerance) const
{
    std::vector<std::unique_ptr<LineString>> outputLines;

    std::unique_ptr<LineString> extRing = conditionRing(poly.getExteriorRing(), tolerance);
    outputLines.emplace_back(extRing.release());

    for (std::size_t i = 0; i < poly.getNumInteriorRing(); i++) {
        auto intRing = conditionRing(poly.getInteriorRingN(i), tolerance);
        outputLines.emplace_back(intRing.release());
    }
    return poly.getFactory()->createMultiLineString(std::move(outputLines));
}


/* private */
std::unique_ptr<LineString>
CoordinateConditioner::conditionRing(
    const LinearRing* ring,
    double tolerance) const
{
    // Convert to working coordinate sequence
    std::unique_ptr<CoordinateSequence> coords = ring->getCoordinatesRO()->clone();
    // Remove duplicate end point from ring
    coords->pop_back();

    // Remove spikes, zigzags and very short segments
    std::unique_ptr<CoordinateSequence> cleaned = removeSmallErrors(*coords, tolerance);

    // Reorient the ring if the start/end is
    // acute. Avoids some stair-stepping effects in interior
    // paths.
    if (startPointAngle(*cleaned) < M_PI / 2.0) {
        std::size_t widestIndex = findWidestAngle(*coords);
        cleaned = reorientCoordinates(*cleaned, widestIndex);
    }

    // Return conditioned LineString
    return ring->getFactory()->createLineString(std::move(cleaned));
}


/* private */
std::unique_ptr<CoordinateSequence>
CoordinateConditioner::removeSmallErrors(
    const CoordinateSequence& coords,
    double tolerance) const
{
    std::size_t sz = coords.size();
    std::unique_ptr<CoordinateSequence> vcs(new CoordinateSequence(0, false, false));
    // vcs->reserve(sz);

    std::size_t i = 0;
    while (i < sz) {

        std::size_t i0 = (i+0) % sz;
        std::size_t i1 = (i+1) % sz;
        std::size_t i2 = (i+2) % sz;
        std::size_t i3 = (i+3) % sz;

        const CoordinateXY& p0 = coords.getAt<CoordinateXY>(i0);
        const CoordinateXY& p1 = coords.getAt<CoordinateXY>(i1);
        const CoordinateXY& p2 = coords.getAt<CoordinateXY>(i2);
        const CoordinateXY& p3 = coords.getAt<CoordinateXY>(i3);

        vcs->add(p0, true);

        if (isZigZag(p0, p1, p2, p3, tolerance)) {
            i += 4;
        }
        else if (isSpike(p0, p1, p2, tolerance)) {
            i += 3;
        }
        else if (isRepeated(p0, p1, tolerance)) {
            i += 2;
        }
        else {
            i += 1;
        }
    }

    return vcs;
}


/* private */
bool
CoordinateConditioner::isZigZag(
    const CoordinateXY& p0,
    const CoordinateXY& p1,
    const CoordinateXY& p2,
    const CoordinateXY& p3,
    double tolerance) const
{
    const double SpikeAngle = M_PI / 360.0;
    double angle0 = Angle::angleBetween(p0, p1, p2);
    double angle1 = Angle::angleBetween(p1, p2, p3);
    double distance = p0.distance(p3);
    if (angle0 <= SpikeAngle && angle1 <= SpikeAngle && distance <= tolerance)
        return true;
    else
        return false;
}


/* private */
bool
CoordinateConditioner::isSpike(
    const CoordinateXY& p0,
    const CoordinateXY& p1,
    const CoordinateXY& p2,
    double tolerance) const
{
    const double SpikeAngle = M_PI / 360.0;
    double angle = Angle::angleBetween(p0, p1, p2);
    double distance = p0.distance(p2);
    if (angle <= SpikeAngle && distance <= tolerance)
        return true;
    else
        return false;
}


/* private */
bool
CoordinateConditioner::isRepeated(
    const CoordinateXY& p0,
    const CoordinateXY& p1,
    double tolerance) const
{
    return p0.distance(p1) <= tolerance;
}


/* private */
std::size_t
CoordinateConditioner::findWidestAngle(
    const CoordinateSequence& coords) const
{
    std::size_t widest = 0;
    double widestAngle = 0.0;
    std::size_t sz = coords.size();
    for (std::size_t i = 0; i < sz; i++) {
        const CoordinateXY &p0 = coords.getAt<CoordinateXY>((i+0) % sz);
        const CoordinateXY &p1 = coords.getAt<CoordinateXY>((i+1) % sz);
        const CoordinateXY &p2 = coords.getAt<CoordinateXY>((i+2) % sz);
        double angle = Angle::angleBetween(p0, p1, p2);
        if (angle > widestAngle) {
            widestAngle = angle;
            widest = i;
        }
    }
    return widest;
}


std::unique_ptr<CoordinateSequence>
CoordinateConditioner::reorientCoordinates(
    const CoordinateSequence& origCoords,
    std::size_t startIndex) const
{
    std::unique_ptr<CoordinateSequence> newCoords(new CoordinateSequence(0, false, false));

    // Start at widest and take everything up to the end
    for (std::size_t i = startIndex; i < origCoords.size(); i++) {
        newCoords->add(origCoords.getAt<CoordinateXY>(i), true);
    }
    // Tack on the rest
    for (std::size_t i = 0; i < startIndex; i++) {
        newCoords->add(origCoords.getAt<CoordinateXY>(i), true);
    }

    return newCoords;
}


double
CoordinateConditioner::startPointAngle(
    const CoordinateSequence& cs) const
{
    const CoordinateXY& p0 = cs.getAt(cs.size()-1);
    const CoordinateXY& p1 = cs.getAt(0);
    const CoordinateXY& p2 = cs.getAt(1);
    return Angle::angleBetween(p0, p1, p2);
}


CoordinateXY
CoordinateConditioner::pointAlong(
    const CoordinateXY& p0,
    const CoordinateXY& p1,
    double segmentLengthFraction)
{
    return CoordinateXY(
        p0.x + segmentLengthFraction * (p1.x - p0.x),
        p0.y + segmentLengthFraction * (p1.y - p0.y));
}


CoordinateXY
CoordinateConditioner::pointAlongDistance(
    const CoordinateXY& p0,
    const CoordinateXY& p1,
    double distance)
{
    double fraction = distance / p0.distance(p1);
    return pointAlong(p0, p1, fraction);
}


/* public static */
std::unique_ptr<MultiLineString>
CoordinateConditioner::densify(const MultiLineString& mls,
    double maxLen)
{
    CoordinateConditioner cc;
    return cc.densifyMultiLineString(mls, maxLen);
}


/* private */
std::unique_ptr<MultiLineString>
CoordinateConditioner::densifyMultiLineString(
    const MultiLineString& mls,
    double maxLen) const
{
    std::vector<std::unique_ptr<LineString>> outputLines;
    for (const auto& lsp : mls) {
        auto ls = static_cast<const LineString*>(lsp.get());
        const CoordinateSequence* coords = ls->getCoordinatesRO();
        std::unique_ptr<CoordinateSequence> densified = densifyCoordinateSequence(coords, maxLen);
        std::unique_ptr<LineString> ols = mls.getFactory()->createLineString(std::move(densified));
        outputLines.emplace_back(ols.release());
    }
    return mls.getFactory()->createMultiLineString(std::move(outputLines));
}

/* Private */
std::unique_ptr<CoordinateSequence>
CoordinateConditioner::densifyCoordinateSequence(
    const CoordinateSequence* coords,
    double maxLen) const
{
    std::size_t sz = coords->size();
    std::unique_ptr<CoordinateSequence> denseCoords(new CoordinateSequence(0, false, false));

    double remainder = 0.0;
    for (std::size_t i = 0; i < sz; i++) {
        auto p0 = coords->getAt<CoordinateXY>((i+0) % sz);
        auto p1 = coords->getAt<CoordinateXY>((i+1) % sz);
        double segLen = p0.distance(p1);

        // First coordinate is always added
        denseCoords->add(p0, true);

        // No room to densify this small segment
        if (maxLen >= segLen)
            continue;

        // Add first sub-point
        denseCoords->add(pointAlong(p0, p1, remainder / segLen), true);
        remainder = segLen - remainder;

        // Add as many other points as we can
        while (remainder > maxLen) {
            remainder -= maxLen;
            denseCoords->add(pointAlong(p0, p1, 1.0 - (remainder/segLen)), true);
        }
    }

    return denseCoords;
}



} // namespace geos.operation.skeletonize
} // namespace geos.operation
} // namespace geos
        
