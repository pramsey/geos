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
 class GeometryFactory;
 class LineString;
 class LineSegment;
 class MultiLineString;
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

// class GEOS_DLL SegmentGraph {
class SegmentGraph {

    using CoordinateXY = geos::geom::CoordinateXY;
    using Geometry = geos::geom::Geometry;
    using GeometryFactory = geos::geom::GeometryFactory;
    using MultiLineString = geos::geom::MultiLineString;
    using LineString = geos::geom::LineString;
    using LineSegment = geos::geom::LineSegment;
    using GeometryLocation = geos::operation::distance::GeometryLocation;

public:

    SegmentGraph(
        std::vector<const Geometry*>& segs,
        std::vector<GeometryLocation>& inoutsegs,
        const GeometryFactory* factory)
        : m_containedEdges(segs)
        , m_inoutEdges(inoutsegs)
        , m_factory(factory)
        {};

    friend std::ostream& operator<< (std::ostream& os, const SegmentGraph& ss);

    std::ostream& toString(std::ostream& os) const;

    // Reads the inputSegments and builds the vertex
    // map/list/adjacency structures
    void build();

    // Finds the longest pairwise path in the graph
    std::unique_ptr<MultiLineString> longestPathSkeleton();

    std::unique_ptr<MultiLineString> shortestPathSkeleton();



    // Returns ordered vector of vertex numbers
    // of the shortest path between start and end vertices
    std::pair<std::vector<uint32_t>, double> shortestPath(uint32_t startVertex, uint32_t endVertex);

    // Returns the list of vertices with cardinality of one
    std::vector<uint32_t> endVertices();

    // Calculate all paths from start vertex and return
    // the one with the highest cost
    std::pair<std::vector<uint32_t>, double> longestPath(uint32_t startVertex, std::vector<uint32_t>& ends);

    // Calculate the overall longest path between all
    // possible pairs of path ends
    std::pair<std::vector<uint32_t>, double> longestPath(std::vector<uint32_t>& ends);

private:


    // All the voronoi edges that are fully contained
    // in the polygon are in this list
    std::vector<const Geometry*>& m_containedEdges;

    // The closest voronoi edge to each in/out point
    // has been noted in this list
    std::vector<GeometryLocation>& m_inoutEdges;

    // Factory to build output geometry
    const GeometryFactory* m_factory;

    // A list of vertices, with one vertex number
    // for each input/output point we matched up
    std::vector<uint32_t> m_inoutVertexList;

    // A list of edges, with the vertex number of each
    // end, and the length of the edge.
    std::vector<std::tuple<uint32_t, uint32_t, double>> m_edgeList;

    // Map (Vertex, Edge)
    // Lookup an edge entry from a vertex number
    std::map<uint32_t, uint32_t> m_edgeMap;

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

    // True if the given coordinate has already
    // been added to the vertex map
    bool isVertex(const CoordinateXY& v) const;

    // Utility functions to break up processing
    // flow in build()
    void buildContainedEdgeList();
    void buildInOutEdgeList();
    void buildAdjacencyList();

    std::unique_ptr<LineString> pathToGeometry(
        std::vector<uint32_t>& vertexPath) const;

    std::unique_ptr<LineString> pathToGeometry(
        uint32_t v0, uint32_t v1) const;

    std::unique_ptr<MultiLineString> mergePairSet(
        std::set<std::pair<uint32_t, uint32_t>>& pairs);

};




} // namespace geos.operation.skeletonize
} // namespace geos.operation
} // namespace geos
        
