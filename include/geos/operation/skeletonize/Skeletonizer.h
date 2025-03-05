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
class GeometryFactory;
class LinearRing;
class LineString;
class MultiLineString;
class MultiPoint;
class Polygon;
}
namespace operation {
namespace distance {
class GeometryLocation;
}
}
}

namespace geos {
namespace operation { // geos::operation
namespace skeletonize { // geos::operation::skeletonize

class GEOS_DLL Skeletonizer {

    using CoordinateSequence = geos::geom::CoordinateSequence;
    using Geometry = geos::geom::Geometry;
    using GeometryFactory = geos::geom::GeometryFactory;
    using MultiLineString = geos::geom::MultiLineString;
    using MultiPoint = geos::geom::MultiPoint;
    using Polygon = geos::geom::Polygon;
    using LinearRing = geos::geom::LinearRing;
    using LineString = geos::geom::LineString;
    using GeometryLocation = geos::operation::distance::GeometryLocation;

public:

    Skeletonizer(const Polygon &poly, const MultiPoint *pts);
    Skeletonizer(const Polygon &poly);

    static std::unique_ptr<MultiLineString> skeletonize(
        const Polygon &poly,
        const MultiPoint &pts);

    static std::unique_ptr<MultiLineString> skeletonize(
        const Polygon& poly);

    std::unique_ptr<MultiLineString> skeletonize();


private:

    const geom::Polygon& inputPolygon;
    const geom::MultiPoint* inputPoints = nullptr;
    const geom::GeometryFactory* inputFactory;
    // double inputWidth;
    // double inputHeight;

    struct SegmentStatistics {
        double numSegments = 0.0;
        double averageLength = 0.0;
        double stdevLength = 0.0;
        double minLength = 0.0;
        double maxLength = 0.0;
        double width = 0.0;
        double height = 0.0;
    };

    void calculateStatistics(
        const CoordinateSequence* cs,
        SegmentStatistics& stats) const;

    void calculateStatistics(
        const Polygon& poly,
        SegmentStatistics& stats) const;

    bool hasInOutPoints() const;

    std::vector<const Geometry*> findContainedEdges(
        const MultiLineString& allEdges) const;

    std::vector<GeometryLocation> findInputOutputEdges(
        const MultiLineString& allEdges) const;

    std::unique_ptr<MultiLineString>
        getGeometry(std::vector<GeometryLocation>& inoutEdges) const;

    std::unique_ptr<MultiLineString>
        getGeometry(std::vector<const Geometry*>& eedges) const;
};


} // namespace geos::operation::skeletonize
} // namespace geos::operation
} // namespace geos
