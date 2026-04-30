/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2026 Paul Ramsey
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#pragma once

#include <geos/export.h>
#include <geos/geom/Coordinate.h>
#include <vector>

// Forward declarations
namespace geos {
namespace geom {
class Curve;
}
}

namespace geos {
namespace operation { // geos::operation
namespace shortest { // geos::operation::shortest

/** \brief
 * Computes shortest paths in a network of Curves using Dijkstra's algorithm.
 *
 * The input is a vector of Curves (LineString, CircularString, CompoundCurve)
 * which form the edges of a graph. Edge weights are curve lengths.
 * The nodes of the graph are the endpoints of the Curves.
 * Coincident endpoints are treated as the same node.
 *
 * Two modes are supported:
 * - Simple: find the shortest path between two supplied points
 *   (snapped to nearest nodes).
 * - Longest shortest path: find the "diameter" path by a triple-Dijkstra sweep.
 */
class GEOS_DLL ShortestPath {

public:

    /** \brief
     * Computes the shortest path between two points in the curve network.
     *
     * The start and end points are snapped to the nearest graph node
     * (curve endpoint) before running Dijkstra's algorithm.
     *
     * @param curves Input vector of Curves forming the network.
     * @param start  Start coordinate (snapped to nearest node).
     * @param end    End coordinate (snapped to nearest node).
     * @param result Output vector of size_t, same size as curves.
     *               0 if the edge is not on the shortest path.
     *               n > 0 gives the position of the edge in the path
     *               (1 = first edge leaving start, 2 = second, etc.).
     */
    static void shortestPath(
        const std::vector<const geom::Curve*>& curves,
        const geom::CoordinateXY& start,
        const geom::CoordinateXY& end,
        std::vector<std::size_t>& result);

    /** \brief
     * Computes the longest shortest path ("diameter") of the curve network.
     *
     * Uses a triple-Dijkstra sweep:
     * 1. Find node A furthest from the centroid of all nodes.
     * 2. Find node B furthest from A.
     * 3. Find node C furthest from B.
     * 4. Return the shortest path from A to C.
     *
     * @param curves Input vector of Curves forming the network.
     * @param result Output vector of size_t, same size as curves.
     *               0 if the edge is not on the path.
     *               n > 0 gives the position of the edge in the path.
     */
    static void longestShortestPath(
        const std::vector<const geom::Curve*>& curves,
        std::vector<std::size_t>& result);

};

} // namespace geos::operation::shortest
} // namespace geos::operation
} // namespace geos
