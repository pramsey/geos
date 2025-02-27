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

 #pragma once

 #include <geos/geom/Coordinate.h>
 // #include <geos/export.h>
 
 #include <memory>
 #include <vector>
 
 // Forward declarations
 namespace geos {
 namespace geom {
 class Geometry;
 class CoordinateSequence;
 class Polygon;
 class LinearRing;
 class LineString;
 class MultiLineString;
 }
 }

namespace geos {
namespace operation {   // geos.operation
namespace skeletonize { // geos.operation.skeletonize

// class GEOS_DLL SegmentGraph {
class CoordinateConditioner {

    using CoordinateXY = geos::geom::CoordinateXY;
    using Geometry = geos::geom::Geometry;
    using CoordinateSequence = geos::geom::CoordinateSequence;
    using Polygon = geos::geom::Polygon;
    using LinearRing = geos::geom::LinearRing;
    using LineString = geos::geom::LineString;
    using MultiLineString = geos::geom::MultiLineString;

public:

    CoordinateConditioner() {};

    static std::unique_ptr<MultiLineString> condition(
        const Polygon& poly,
        double tolerance);

    static std::unique_ptr<MultiLineString> densify(
        const MultiLineString& mls,
        double maxLen);

    static CoordinateXY pointAlong(
        const CoordinateXY& p0,
        const CoordinateXY& p1,
        double segmentLengthFraction);

    static CoordinateXY pointAlongDistance(
        const CoordinateXY& p0,
        const CoordinateXY& p1,
        double distance);


private:

    std::unique_ptr<MultiLineString> conditionPolygon(
        const Polygon& poly,
        double tolerance) const;

    std::unique_ptr<LineString> conditionRing(
        const LinearRing* lr,
        double tolerance) const;

    std::unique_ptr<CoordinateSequence> removeSmallErrors(
        const CoordinateSequence& coords,
        double tolerance) const;

    bool isRepeated(
        const CoordinateXY& p0,
        const CoordinateXY& p1,
        double tolerance) const;

    bool isSpike(
        const CoordinateXY& p0,
        const CoordinateXY& p1,
        const CoordinateXY& p2,
        double tolerance) const;

    bool isZigZag(
        const CoordinateXY& p0,
        const CoordinateXY& p1,
        const CoordinateXY& p2,
        const CoordinateXY& p3,
        double tolerance) const;

    std::size_t findWidestAngle(
        const CoordinateSequence& coords) const;

    std::unique_ptr<CoordinateSequence> reorientCoordinates(
        const CoordinateSequence& origCoords,
        std::size_t startIndex) const;

    double startPointAngle(
        const CoordinateSequence& cs) const;

    std::unique_ptr<MultiLineString> densifyMultiLineString(
        const MultiLineString& mls,
        double maxLen) const;

    std::unique_ptr<CoordinateSequence> densifyCoordinateSequence(
        const CoordinateSequence* cs,
        double maxLen) const;

};




} // namespace geos.operation.skeletonize
} // namespace geos.operation
} // namespace geos
        
