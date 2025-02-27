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
 #include <geos/geom/GeometryFactory.h>
 #include <geos/geom/LineSegment.h>
 #include <geos/geom/LineString.h>
 #include <geos/geom/MultiLineString.h>
 #include <geos/util/GEOSException.h>
 
 #include <queue>


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
SegmentGraph::build(void)
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

void
SegmentGraph::clear(void)
{
    m_adj.clear();
    m_vertexList.clear();
    m_vertexMap.clear();
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


/* public */
std::vector<uint32_t>
SegmentGraph::endVertices()
{
    std::vector<uint32_t> ends;
    // Map (Coordinate, (Id, Cardinality))
    for (const auto& entry : m_vertexMap) {
        uint32_t id = entry.second.first;
        uint32_t cardinality = entry.second.second;
        if (cardinality == 1) {
            ends.push_back(id);
        }
    }
    return ends;
}


std::pair<std::vector<uint32_t>, double>
SegmentGraph::shortestPath(uint32_t startVertex, uint32_t endVertex)
{
    const double UnknownDist = std::numeric_limits<double>::max();
    const uint32_t UnknownVertex = std::numeric_limits<uint32_t>::max();
    std::vector<double> dist(m_vertexCount, UnknownDist);
    std::vector<uint32_t> parent(m_vertexCount, UnknownVertex); // To reconstruct path

    std::priority_queue<
        std::pair<double, uint32_t>,
        std::vector<std::pair<double, uint32_t>>,
        std::greater<>> pq;

    dist[startVertex] = 0.0;
    pq.emplace(0.0, startVertex);

    while (!pq.empty()) {
        double d = pq.top().first;
        uint32_t u = pq.top().second;
        pq.pop();

        // Stop early if endVertex is reached
        if (u == endVertex) break;

        if (d > dist[u]) continue;

        for (const auto& pr : m_adj[u]) {
            uint32_t v = pr.first;
            double weight = pr.second;
            double newDist = dist[u] + weight;
            if (newDist < dist[v]) {
                dist[v] = newDist;
                parent[v] = u; // Store where v came from
                pq.emplace(newDist, v);
            }
        }
    }

    // If target is unreachable
    if (std::isnan(dist[endVertex]))
        return {{}, UnknownDist};

    // Reconstruct path from target to source
    std::vector<uint32_t> path;
    for (uint32_t at = endVertex; at != UnknownVertex; at = parent[at]) {
        path.push_back(at);
    }
    reverse(path.begin(), path.end());

    return {path, dist[endVertex]};
}


std::vector<uint32_t>
SegmentGraph::longestPath(uint32_t startVertex, std::vector<uint32_t>& ends)
{
    double maxDist = 0.0;
    std::vector<uint32_t> maxPath;
    for (uint32_t endVertex : ends) {
        if (endVertex == startVertex)
            continue;
        auto result = shortestPath(startVertex, endVertex);
        double dist = result.second;
        if (dist > maxDist) {
            maxDist = dist;
            maxPath = result.first;
        }
    }
    return maxPath;
}


std::vector<uint32_t>
SegmentGraph::longestPath(std::vector<uint32_t>& ends)
{
    // Pick any vertex in the set of ends and use that
    // to calculate a temporary maximal route
    uint32_t aStart = ends[0];
    std::vector<uint32_t> longishPath = longestPath(aStart, ends);
    // Taking the other end of the maximal route, and finding
    // the longest route from there, gives us the
    // overall longest route for the set
    uint32_t longerStart = longishPath[longishPath.size()-1];
    std::vector<uint32_t> longPath = longestPath(longerStart, ends);
    return longPath;
}


std::unique_ptr<LineString>
SegmentGraph::longestPath()
{
    // Get fresh build of graph structures
    clear();
    build();

    std::vector<uint32_t> graphEnds = endVertices();
    std::vector<uint32_t> vertexPath = longestPath(graphEnds);

    // Convert path of ids into geometric path
    std::unique_ptr<CoordinateSequence> cs = std::make_unique<CoordinateSequence>();
    for (auto vertexId : vertexPath) {
        CoordinateXY c = m_vertexList[vertexId];
        cs->add(c, false);
    }
    return m_inputSegments.getFactory()->createLineString(std::move(cs));
}


std::unique_ptr<MultiLineString>
SegmentGraph::longestPaths()
{
    // Get fresh build of graph structures
    clear();
    build();

    std::vector<uint32_t> graphEnds = endVertices();
    std::vector<uint32_t> vertexPath = longestPath(graphEnds);

    // Convert path of ids into geometric path
    std::unique_ptr<CoordinateSequence> cs = std::make_unique<CoordinateSequence>();
    for (auto vertexId : vertexPath) {
        CoordinateXY c = m_vertexList[vertexId];
        cs->add(c, false);
    }
    auto path = m_inputSegments.getFactory()->createLineString(std::move(cs));
    std::vector<std::unique_ptr<LineString>> lines;
    lines.emplace_back(path.release());
    return m_inputSegments.getFactory()->createMultiLineString(std::move(lines));
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
        
