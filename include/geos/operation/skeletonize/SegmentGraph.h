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
    void buildAdjacencyList();


private:

    const MultiLineString& m_inputSegments;

    // Map (Vertex, (ID, Cardinality))
    // Iterate to find all dangling vertices
    // Lookup to find given version identifier
    std::map<CoordinateXY, std::pair<uint32_t, uint32_t>> m_vertexMap;

    std::vector<CoordinateXY> m_vertexList;

    // Vector (ID0, (ID1, Weight)[])
    // Vertex adjacency list, m_vertexCount in length, each
    // entry with a list of adjacent vertices and cost to
    // that vertex
    std::vector<std::vector<std::pair<uint32_t, double>>> m_adj;

    // Number of unique vertices in the collection
    uint32_t m_vertexCount = 0;



    uint32_t mapVertex(const CoordinateXY& v);

};




} // namespace geos.operation.skeletonize
} // namespace geos.operation
} // namespace geos
        
