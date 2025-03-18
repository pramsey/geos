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
class LineString;
class MultiLineString;
class MultiPoint;
class Point;
}
namespace operation {
namespace distance {
class GeometryLocation;
}
}
}

namespace geos {
namespace operation {   // geos.operation
namespace skeletonize { // geos.operation.skeletonize

typedef std::map<const geos::geom::Geometry*, std::vector<geos::operation::distance::GeometryLocation>> InputOutputMap;

class InputOutputs {

    using CoordinateXY = geos::geom::CoordinateXY;
    using Geometry = geos::geom::Geometry;
    using CoordinateSequence = geos::geom::CoordinateSequence;
    using LineString = geos::geom::LineString;
    using MultiLineString = geos::geom::MultiLineString;
    using MultiPoint = geos::geom::MultiPoint;
    using Point = geos::geom::Point;
    using GeometryLocation = geos::operation::distance::GeometryLocation;

public:

    InputOutputs() {};

    static std::unique_ptr<MultiLineString> addInputOutputGaps(
        const MultiLineString& mls,
        const std::vector<const Point*>& pts,
        double tolerance);


private:

    std::unique_ptr<MultiLineString> process(
        const MultiLineString& mls,
        const std::vector<const Point*>& pts,
        double tolerance) const;

    std::size_t coordIndex(
        const GeometryLocation& loc) const;

    const CoordinateSequence* getCoordinates(
        const GeometryLocation& loc) const;

    void addGappedPair(
        const GeometryLocation& loc,
        double tolerance,
        CoordinateSequence& newCoords) const;

    InputOutputMap findClosestLocations(
        const MultiLineString& mls,
        const std::vector<const Point*>& pts) const;
};




} // namespace geos.operation.skeletonize
} // namespace geos.operation
} // namespace geos
        
