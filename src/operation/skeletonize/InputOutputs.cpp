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

#include <geos/operation/skeletonize/InputOutputs.h>
#include <geos/operation/skeletonize/CoordinateConditioner.h>
#include <geos/operation/distance/GeometryLocation.h>
#include <geos/algorithm/Angle.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/LineString.h>
#include <geos/geom/MultiLineString.h>
#include <geos/geom/MultiPoint.h>
#include <geos/geom/Point.h>
#include <geos/operation/distance/IndexedFacetDistance.h>

#include <queue>


using geos::algorithm::Angle;
using geos::geom::CoordinateXY;
using geos::geom::CoordinateSequence;
using geos::geom::Geometry;
using geos::geom::GeometryFactory;
using geos::geom::LineString;
using geos::geom::MultiLineString;
using geos::geom::MultiPoint;
using geos::geom::Point;
using geos::operation::distance::GeometryLocation;
using geos::operation::distance::IndexedFacetDistance;


namespace geos {
namespace operation {   // geos.operation
namespace skeletonize { // geos.operation.skeletonize


/* public static */
std::unique_ptr<MultiLineString>
InputOutputs::addInputOutputGaps(
    const MultiLineString& mls,
    const MultiPoint& pts,
    double tolerance)
{
    InputOutputs iop;
    return iop.process(mls, pts, tolerance);
}

/*
 *
 * Convert output of conditioning function to MultiLineString
 * Build IndexedFacetDistance on MultiLineString
 * For each MultiPoint input point, find GeometryLocation on MLS
 * add each GL to a map<CoordinateSeq*, vector<SegNum>>
 * for each entry in the map, build a new coordinateSeq with the
 * special point separation
 *
 */

const CoordinateSequence*
InputOutputs::getCoordinates(
    const GeometryLocation& loc) const
{
    const Geometry* geomComp = loc.getGeometryComponent();
    const LineString* ls = static_cast<const LineString*>(geomComp);
    return ls->getCoordinatesRO();
}


std::size_t
InputOutputs::coordIndex(
    const GeometryLocation& loc) const
{
    const CoordinateSequence* coords = getCoordinates(loc);
    std::size_t segmentIndex = loc.getSegmentIndex();
    auto c0 = coords->getAt<CoordinateXY>(segmentIndex);
    auto c1 = coords->getAt<CoordinateXY>(segmentIndex+1);
    auto p0 = loc.getCoordinate();
    double d0 = c0.distance(p0);
    double d1 = c1.distance(p0);
    return d0 < d1 ? segmentIndex : segmentIndex+1;
}


void
InputOutputs::addGappedPair(
    const GeometryLocation& loc,
    double tolerance,
    CoordinateSequence& newCoords) const
{
    const CoordinateSequence* coords = getCoordinates(loc);
    std::size_t segmentIndex = loc.getSegmentIndex();
    std::size_t sz = coords->size();
    std::size_t s0 = (segmentIndex + (sz-1)) % sz;
    std::size_t s2 = (segmentIndex + (sz+1)) % sz;

    auto p0 = coords->getAt<CoordinateXY>(s0);
    auto p1 = loc.getCoordinate();
    auto p2 = coords->getAt<CoordinateXY>(s2);

    auto n0 = CoordinateConditioner::pointAlongDistance(p1, p0, tolerance/10);
    auto n2 = CoordinateConditioner::pointAlongDistance(p1, p2, tolerance/10);

    newCoords.add(n0, true);
    newCoords.add(n2, true);
    return;
}


/* private */
InputOutputMap
InputOutputs::findClosestLocations(
    const MultiLineString& mls,
    const MultiPoint& pts) const
{
    IndexedFacetDistance ifd(&mls);
    InputOutputMap ioMap;

    //
    // For each point of interest, find the nearest point on the
    // target LineString and store the GeometryLocation in the map<>.
    //
    // std::cout << "InputOutputs::process " << std::endl;
    // std::cout << " iterating on components of MultiPoint" << std::endl;
    for (std::size_t i = 0; i < pts.getNumGeometries(); i++) {
        const Point* pt = pts.getGeometryN(i);
        std::vector<GeometryLocation> locs = ifd.nearestLocations(pt);
        // std::cout << "  i = " << i << std::endl;
        // std::cout << "    locs[0] = " << locs[0].toString() << std::endl;
        // std::cout << "    locs[1] = " << locs[1].toString() << std::endl;

        // Add the found location to the map associated with the component it was found in
        GeometryLocation& gl = locs[0];
        std::size_t ptIndex = coordIndex(gl);
        const Geometry* geomComp = gl.getGeometryComponent();
        auto ptCoord = pt->getCoordinatesRO()->getAt<CoordinateXY>(0);
        ioMap[geomComp].emplace_back(geomComp, ptIndex, ptCoord);
    }

    return ioMap;
}


/* private */
std::unique_ptr<MultiLineString>
InputOutputs::process(
    const MultiLineString& mls,
    const MultiPoint& pts,
    double tolerance) const
{
    (void)tolerance; // xxxxxx

    //
    // For each input/output point, find the segment of the
    // multilinestring that is closest, so that later we can
    // replace that segment with a gapped pair of points
    // that will drive the voronoi diagram to put an edge
    // very close to the inout/output point
    //
    auto ioMap = findClosestLocations(mls, pts);

    //
    // For each LineString in the input MultiLineString, iterate and
    // replace the targetted points with gapped points around the
    // input/output points.
    //
    // std::cout << std::endl;
    // std::cout << " iterating on components of MultiLineString" << std::endl;
    std::vector<std::unique_ptr<LineString>> outputLines;
    for (std::size_t i = 0; i < mls.getNumGeometries(); i++) {

        const LineString* ls = mls.getGeometryN(i);
        const CoordinateSequence* coords = ls->getCoordinatesRO();

        //
        // No input/output points on this ring.
        // Take a copy and move on.
        //
        if (ioMap.find(ls) == ioMap.end()) {
            outputLines.emplace_back(ls->clone().release());
            continue;
        }

        std::vector<GeometryLocation>& geomLocs = ioMap[ls];

        //
        // We want to process the locations in order so we can
        // fill in the coordinates we are retaining. Pre-sort
        // on index number.
        //
        std::sort(geomLocs.begin(), geomLocs.end(),
            [](const GeometryLocation &a, const GeometryLocation &b) {
                return a.getSegmentIndex() < b.getSegmentIndex();
            });

        std::unique_ptr<CoordinateSequence> newSeq(new CoordinateSequence(0, false, false));

        std::size_t start = 0;
        for (const GeometryLocation& loc : geomLocs) {
            // Fill up to the next replacement coordinate
            std::size_t idx = loc.getSegmentIndex();
            for (std::size_t j = start; j < idx; j++) {
                newSeq->add(coords->getAt<CoordinateXY>(j), true);
            }
            // Replace it with a gapped pair
            addGappedPair(loc, tolerance, *newSeq);
            start = idx+1;
        }
        for (std::size_t j = start; j < coords->size(); j++) {
            newSeq->add(coords->getAt<CoordinateXY>(j), true);
        }

        // std::cout << "input coords" << std::endl;
        // std::cout << "LINESTRING" << *coords << std::endl;
        // std::cout << "gapped CoordinateSequence" << std::endl;
        // std::cout << "LINESTRING" << *newSeq << std::endl;

        auto newLs = mls.getFactory()->createLineString(std::move(newSeq));
        outputLines.emplace_back(newLs.release());
    }

    return mls.getFactory()->createMultiLineString(std::move(outputLines));
}


} // namespace geos.operation.skeletonize
} // namespace geos.operation
} // namespace geos
        
