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

 #include <geos/operation/skeletonize/SegmentGraph.h>
 #include <geos/geom/Coordinate.h>
 #include <geos/geom/CoordinateSequence.h>
 #include <geos/geom/LineSegment.h>
 #include <geos/geom/LineString.h>
 #include <geos/geom/MultiLineString.h>
 #include <geos/util/GEOSException.h>
 
 
 using geos::geom::Coordinate;
 using geos::geom::CoordinateSequence;
 using geos::geom::LineSegment;
 using geos::geom::LineString;
 using geos::geom::MultiLineString;
 

namespace geos {
namespace operation {   // geos.operation
namespace skeletonize { // geos.operation.skeletonize


/* private */
void
SegmentGraph::buildAdjacencyList(void)
{
    std::vector<std::tuple<uint32_t, uint32_t, double>> edges;

    // The inputSegments are the MultiLineString output
    // of VoronoiDiagramBuilder and are expected to all
    // be two-point LineString with shared end points.
    for (uint32_t i = 0; i < m_inputSegments.getNumGeometries(); i++) {
        const LineString *ls = m_inputSegments.getGeometryN(i);
        const CoordinateSequence *cs = ls->getCoordinatesRO();
        auto& p0 = cs->getAt<CoordinateXY>(0);
        auto& p1 = cs->getAt<CoordinateXY>(cs->size()-1);

        // Add start/end indexes to a temporary list
        uint32_t i0 = mapVertex(p0);
        uint32_t i1 = mapVertex(p1);
        double weight = p0.distance(p1);
        edges.emplace_back(i0, i1, weight);
    }

    // Populate adjacency list
    m_adj.resize(m_vertexCount);
    for (const auto& edge : edges) {
        uint32_t i0 = std::get<0>(edge);
        uint32_t i1 = std::get<1>(edge);
        double weight = std::get<2>(edge);
        m_adj[i0].emplace_back(i1, weight);
        m_adj[i1].emplace_back(i0, weight);
    }

    // Populate a list of vertex locations so we can
    // build up spatial linestrings from sets of vertex
    // path ids after running the algorithm.
    m_vertexList.resize(m_vertexCount);
    for (const auto& entry : m_vertexMap) {
        uint32_t id = entry.second.first;
        const CoordinateXY& v = entry.first;
        m_vertexList[id] = v;
    }

    return;
}


/* private */
uint32_t
SegmentGraph::mapVertex(const CoordinateXY& v)
{
    uint32_t count = 1;
    uint32_t id = m_vertexCount;
    auto search = m_vertexMap.find(v);
    // Vertex is already in the map
    if (search != m_vertexMap.end()) {
        id = search->second.first;
        count = search->second.second + 1;
    }
    // Vertex is new to the map
    else {
        id = m_vertexCount++;
    }
    m_vertexMap[v] = {id, count};
    return id;
}


std::ostream&
operator<< (std::ostream& os, const SegmentGraph& ss)
{
    return ss.toString(os);
}

std::ostream&
SegmentGraph::toString(std::ostream& os) const
{
    os << "SegmentGraph " << std::endl;
    for (uint32_t i = 0; i < m_adj.size(); i++) {
        os << " [" << i << ": ";
        for (uint32_t j = 0; j < m_adj[i].size(); j++) {
            if (j) os << ",";
            os << m_adj[i][j].first << " " << m_adj[i][j].second;
        }
        os << "] " << std::endl;;
    }
    return os;
}



} // namespace geos.operation.skeletonize
} // namespace geos.operation
} // namespace geos
        
