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

#include <vector>
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
class MultiPolygon;
class Polygon;
class Point;
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
    using MultiPolygon = geos::geom::MultiPolygon;
    using Point = geos::geom::Point;
    using Polygon = geos::geom::Polygon;
    using LinearRing = geos::geom::LinearRing;
    using LineString = geos::geom::LineString;
    using GeometryLocation = geos::operation::distance::GeometryLocation;

public:

    Skeletonizer(const Geometry* poly, const MultiPoint* pts);
    Skeletonizer(const Geometry* poly);

    static std::unique_ptr<MultiLineString> skeletonize(
        const Geometry* polys,
        const MultiPoint* pts,
        double tolerance = 0.0,
        double conditioningLength = 0.0);

    static std::unique_ptr<MultiLineString> skeletonize(
        const Geometry* polys,
        double tolerance = 0.0,
        double conditioningLength = 0.0);


    double getTolerance() const {
        return m_tolerance;
    }

    void setTolerance(double new_tolerance) {
        m_tolerance = new_tolerance;
    }

    double getConditioningLength() const {
        return m_conditioningLength;
    }

    void setConditioningLength(double new_conditioningLength) {
        m_conditioningLength = new_conditioningLength;
    }

private:

    const geom::Geometry* m_inputGeometry;
    const geom::MultiPoint* m_inputPoints = nullptr;
    const geom::GeometryFactory* m_inputFactory = nullptr;
    double m_tolerance = 0.0;
    double m_conditioningLength = 0.0;

    struct SegmentStatistics {
        double numSegments = 0.0;
        double averageLength = 0.0;
        double stdevLength = 0.0;
        double minLength = 0.0;
        double maxLength = 0.0;
        double width = 0.0;
        double height = 0.0;
        double M2 = 0.0;
        double maxX = 0.0;
        double maxY = 0.0;
    };

    void calculateStatistics(
        const CoordinateSequence* cs,
        SegmentStatistics& stats) const;

    void initializeStatistics(
        SegmentStatistics& stats) const;

    SegmentStatistics calculateStatistics(
        const Polygon& poly) const;

    double defaultTolerance(
        const SegmentStatistics& stats) const;

    double defaultConditioningLength(
        const SegmentStatistics& stats) const;

    std::vector<const Geometry*> findContainedEdges(
        const MultiLineString& allEdges) const;

    std::vector<GeometryLocation> findInputOutputEdges(
        const MultiLineString& allEdges,
        const std::vector<const Point*>& points) const;

    std::unique_ptr<MultiLineString>
        getGeometry(std::vector<GeometryLocation>& inoutEdges) const;

    std::unique_ptr<MultiLineString>
        getGeometry(std::vector<const Geometry*>& edges) const;

    void skeletonizePolygon(
        const Polygon* poly,
        const std::vector<const Point*>& points,
        std::vector<std::unique_ptr<LineString>>& skelLines);

    std::unique_ptr<MultiLineString> skeletonize();
};


} // namespace geos::operation::skeletonize
} // namespace geos::operation
} // namespace geos
