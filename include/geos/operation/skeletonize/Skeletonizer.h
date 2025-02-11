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

#include <geos/export.h>

#include <memory>

// Forward declarations
namespace geos {
namespace geom {
class CoordinateSequence;
class Geometry;
class LinearRing;
class LineString;
class MultiLineString;
class Polygon;
}
}

namespace geos {
namespace operation { // geos::operation
namespace skeletonize { // geos::operation::skeletonize

class GEOS_DLL Skeletonizer {

    using CoordinateSequence = geos::geom::CoordinateSequence;
    using Geometry = geos::geom::Geometry;
    using MultiLineString = geos::geom::MultiLineString;
    using Polygon = geos::geom::Polygon;
    using LinearRing = geos::geom::LinearRing;
    using LineString = geos::geom::LineString;

public:

    Skeletonizer(const Polygon &poly)
        : inputPolygon(poly)
        {};

    std::unique_ptr<LineString> skeletonize();

    static std::unique_ptr<LineString> skeletonize(const Polygon &poly);

    std::unique_ptr<Geometry> densifyDefault(const Polygon* poly, double tolerance);

    std::unique_ptr<Geometry> densifyUniformly(const Polygon* poly, double tolerance);

    std::unique_ptr<LinearRing> densifyUniformly(const LinearRing* ring, double tolerance);


private:

    struct SegmentStatistics {
        std::size_t numSegments;
        double totalLength;
        double medianLength;
        double averageLength;
        double minLength;
        double maxLength;
    };

    const geom::Polygon& inputPolygon;
    // double inputWidth;
    // double inputHeight;

    double startPointAngle(const CoordinateSequence* cs) const;

    std::size_t findWidestAngle(const CoordinateSequence* ring) const;

    std::unique_ptr<CoordinateSequence> startSequenceAtWidest(const CoordinateSequence* seq);

};


} // namespace geos::operation::skeletonize
} // namespace geos::operation
} // namespace geos
