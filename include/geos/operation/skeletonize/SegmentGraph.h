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

 #include <geos/geom/LineSegment.h>
 #include <geos/geom/Coordinate.h>
 // #include <geos/export.h>
 
 #include <memory>
 #include <vector>
 
 // Forward declarations
 namespace geos {
 namespace geom {
 class Geometry;
 class LineString;
 class LineSegment;
 class MultiLineString;
 }
 }

namespace geos {
namespace operation {   // geos.operation
namespace skeletonize { // geos.operation.skeletonize

// class GEOS_DLL SegmentGraph {
class SegmentGraph {

    using CoordinateXY = geos::geom::CoordinateXY;
    using Geometry = geos::geom::Geometry;
    using MultiLineString = geos::geom::MultiLineString;
    using LineString = geos::geom::LineString;
    using LineSegment = geos::geom::LineSegment;

public:

    SegmentGraph(const MultiLineString& segs)
        : m_inputSegments(segs)
        {};

    friend std::ostream& operator<< (std::ostream& os, const SegmentGraph& ss);

    std::ostream& toString(std::ostream& os) const;

    // Reads the inputSegments and builds the vertex
    // map/list/adjacency structures
    void build();

    // Removes all the data structures supporting the graph
    void clear();

    // Based on the input segments, calculate the
    // longest path through the graph
    std::unique_ptr<LineString> longestPath();

    // Based on the input segments, calculate the
    // longest path through the graph
    std::unique_ptr<MultiLineString> longestPaths();

    // Returns ordered vector of vertex numbers
    // of the shortest path between start and end vertices
    std::pair<std::vector<uint32_t>, double> shortestPath(uint32_t startVertex, uint32_t endVertex);

    // Returns the list of vertices with cardinality of one
    std::vector<uint32_t> endVertices();

    // Calculate all paths from start vertex and return
    // the one with the highest cost
    std::vector<uint32_t> longestPath(uint32_t startVertex, std::vector<uint32_t>& ends);

    // Calculate the overall longest path between all
    // possible pairs of path ends
    std::vector<uint32_t> longestPath(std::vector<uint32_t>& ends);

private:

    const MultiLineString& m_inputSegments;

    // Map (Vertex, (ID, Cardinality))
    // Iterate to find all dangling vertices
    // Lookup to find version identifier for coordinate
    std::map<CoordinateXY, std::pair<uint32_t, uint32_t>> m_vertexMap;

    // Vector (Coordinate)
    // Constant time lookup of the input coordinate for
    // a given unique vertex identifier
    std::vector<CoordinateXY> m_vertexList;

    // Vector (ID0, (ID1, Weight)[])
    // Vertex adjacency list, m_vertexCount in length, each
    // entry with a list of adjacent vertices and cost to
    // that vertex
    std::vector<std::vector<std::pair<uint32_t, double>>> m_adj;

    // Number of unique vertices in the collection
    uint32_t m_vertexCount = 0;

    // Lookup or generate the unique vertex number for this
    // coordinate
    uint32_t mapVertex(const CoordinateXY& v);

};




} // namespace geos.operation.skeletonize
} // namespace geos.operation
} // namespace geos
        
