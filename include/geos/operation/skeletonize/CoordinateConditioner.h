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


public:

    CoordinateConditioner() {};

    static std::unique_ptr<Geometry> condition(
        const Polygon* poly,
        double tolerance,
        double maxLen);

    std::unique_ptr<Geometry> conditionPolygon(
        const Polygon* poly,
        double tolerance,
        double maxLen) const;


private:

    std::unique_ptr<LinearRing> conditionRing(
        const LinearRing* cs,
        double tolerance,
        double maxLen) const;

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

    CoordinateXY pointAlong(
        const CoordinateXY& p0,
        const CoordinateXY& p1,
        double segmentLengthFraction) const;

    std::unique_ptr<CoordinateSequence> densify(
        const CoordinateSequence& coords,
        double maxLen) const;

};




} // namespace geos.operation.skeletonize
} // namespace geos.operation
} // namespace geos
        
