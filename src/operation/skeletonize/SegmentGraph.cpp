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

 #include <geos/operation/distance/GeometryLocation.h>
 #include <geos/operation/skeletonize/SegmentGraph.h>
 #include <geos/geom/Coordinate.h>
 #include <geos/geom/CoordinateSequence.h>
 #include <geos/geom/GeometryFactory.h>
 #include <geos/geom/LineSegment.h>
 #include <geos/geom/LineString.h>
 #include <geos/geom/MultiLineString.h>
 #include <geos/util/GEOSException.h>
 #include <geos/operation/linemerge/LineMerger.h>

 #include <queue>


 using geos::geom::Coordinate;
 using geos::geom::CoordinateSequence;
 using geos::geom::LineSegment;
 using geos::geom::LineString;
 using geos::geom::MultiLineString;
 using geos::operation::distance::GeometryLocation;
 using geos::operation::linemerge::LineMerger;
 

namespace geos {
namespace operation {   // geos.operation
namespace skeletonize { // geos.operation.skeletonize


/* private */
void
SegmentGraph::buildAdjacencyList()
{
    // Populate adjacency list
    m_adj.resize(m_vertexCount);
    for (const auto& edge : m_edgeList) {
        uint32_t i0 = std::get<0>(edge);
        uint32_t i1 = std::get<1>(edge);
        double weight = std::get<2>(edge);
        // All our edges are bidirectional
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
}


/* private */
void
SegmentGraph::buildContainedEdgeList()
{
    // The containedEdges are the MultiLineString output
    // of VoronoiDiagramBuilder and are expected to all
    // be two-point LineStrings with shared end points.
    for (const Geometry* geom : m_containedEdges) {

        const LineString *ls = dynamic_cast<const LineString*>(geom);
        if(!ls)
            continue;

        // Read start/end coordinates from the edge
        const CoordinateSequence *cs = ls->getCoordinatesRO();
        auto& p0 = cs->getAt<CoordinateXY>(0);
        auto& p1 = cs->getAt<CoordinateXY>(cs->size()-1);

        // Look-up/generate the vertex number for each
        // vertex coordinate, and calculate the edge weight
        uint32_t i0 = mapVertex(p0);
        uint32_t i1 = mapVertex(p1);
        double weight = p0.distance(p1);
        m_edgeList.emplace_back(i0, i1, weight);
    }
}


/* private */
void
SegmentGraph::buildInOutEdgeList()
{
    for (const GeometryLocation& loc : m_inoutEdges) {

        const Geometry* geomComp = loc.getGeometryComponent();
        const LineString *ls = dynamic_cast<const LineString*>(geomComp);
        if(!ls)
            continue;

        // Read start/end coordinates from the edge
        const CoordinateSequence *cs = ls->getCoordinatesRO();
        const auto& p0 = cs->getAt<CoordinateXY>(0);
        const auto& p1 = cs->getAt<CoordinateXY>(cs->size()-1);

        // Read input/output coordinate from the location
        const auto& p2 = loc.getCoordinate();

        // One end of the edge should connect to the
        // existing set of interior edges. (Hope so!)
        // We add the connecting vertex, and the known
        // input/output point, thus ensuring the graph
        // we generate connects up to the input/output
        // points.
        bool h0 = isVertex(p0);
        bool h1 = isVertex(p1);
        if (h0 && h1) {
            throw std::runtime_error("edge doubly connects!");
        }
        else if (!(h0||h1)) {
            throw std::runtime_error("edge does not connect!");
        }
        else {
            // Vertex that connects to graph
            uint32_t i0 = h0 ? mapVertex(p0) : mapVertex(p1);
            // Vertex that is input/output point
            uint32_t i2 = mapVertex(p2);
            double weight = p0.distance(p2);
            m_edgeList.emplace_back(i0, i2, weight);
            // Save the vertex id of the in/out point
            // so we can use as start/end point of shortest
            // path calculations later
            m_inoutVertexList.push_back(i2);
        }
    }
}


/* private */
void
SegmentGraph::build(void)
{
    // Reset all the local data structures
    m_adj.clear();
    m_vertexList.clear();
    m_vertexMap.clear();
    m_edgeMap.clear();
    m_edgeList.clear();
    m_inoutVertexList.clear();

    // Add all edges from m_containedEdges
    buildContainedEdgeList();

    // Add all edges from m_inoutEdges
    buildInOutEdgeList();

    // Build m_edgeList into m_adj and m_vertexList
    buildAdjacencyList();

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


/* private */
bool
SegmentGraph::isVertex(const CoordinateXY& v) const
{
    auto search = m_vertexMap.find(v);
    return search != m_vertexMap.end();
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

    // std::cout << "shortestPath(" <<  startVertex << "," << endVertex << ")" << std::endl;

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


std::pair<std::vector<uint32_t>, double>
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
    return {maxPath, maxDist};
}


std::pair<std::vector<uint32_t>, double>
SegmentGraph::longestPath(std::vector<uint32_t>& ends)
{
    // Pick any vertex in the set of ends and use that
    // to calculate a temporary maximal route
    uint32_t aStart = ends[0];
    std::vector<uint32_t> longishPath = longestPath(aStart, ends).first;
    // Taking the other end of the maximal route, and finding
    // the longest route from there, gives us the
    // overall longest route for the set
    uint32_t longerStart = longishPath[longishPath.size()-1];
    return longestPath(longerStart, ends);
}



std::vector<std::unique_ptr<LineString>>
SegmentGraph::longestPathSkeleton()
{
    // Get fresh build of graph structures
    build();

    // List of cardinality-1 vertices in the graph
    std::vector<uint32_t> graphEnds = endVertices();
    // Longest pairwise path between those vertices
    std::vector<uint32_t> vertexPath = longestPath(graphEnds).first;
    // Wrap LineString to output MultiLineString
    std::vector<std::unique_ptr<LineString>> lines;
    lines.emplace_back(pathToGeometry(vertexPath).release());
    return lines;
}




std::vector<std::unique_ptr<LineString>>
SegmentGraph::shortestPathSkeleton()
{
    // Get fresh build of graph structures
    build();

    // In case of a single in/out vertices we build the longest
    // path that terminates at that point.
    if (m_inoutVertexList.size() == 1) {
        std::vector<uint32_t> graphEnds = endVertices();
        uint32_t startVertex = m_inoutVertexList[0];
        std::pair<std::vector<uint32_t>, double> vertexPath = longestPath(startVertex, graphEnds);
        auto ls = pathToGeometry(vertexPath.first);
        std::vector<std::unique_ptr<LineString>> lines;
        lines.emplace_back(ls.release());
        return lines;
    }

    // For each pairing of in/out vertices, generate
    // the shortest path, and add that to our set
    // of path pairs. We should end up with the unique
    // set of segments (as defined by pairs)
    // that form our skeleton.
    std::vector<uint32_t> pathEnds = m_inoutVertexList;
    std::set<std::pair<uint32_t, uint32_t>> pathPairs;
    while (!pathEnds.empty()) {
        uint32_t start = pathEnds.back();
        pathEnds.pop_back();
        for (uint32_t end : pathEnds) {
            // std::pair<std::vector<uint32_t>,double>
            auto sp = shortestPath(start, end).first;
            for (std::size_t i = 1; i < sp.size(); i++) {
                uint32_t a = sp[i-1], b = sp[i], tmp;
                if (a > b) {
                    tmp = a; a = b; b = tmp;
                }
                pathPairs.emplace(a, b);
            }
        }
    }

    return pairsToGeometry(pathPairs);
}


std::unique_ptr<LineString>
SegmentGraph::pathToGeometry(
    std::vector<uint32_t>& vertexPath) const
{
    auto cs = std::make_unique<CoordinateSequence>();
    for (uint32_t vertex : vertexPath) {
        CoordinateXY c = m_vertexList[vertex];
        cs->add(c, false);
    }
    return m_factory->createLineString(std::move(cs));
}

std::unique_ptr<LineString>
SegmentGraph::pathToGeometry(
    uint32_t v0, uint32_t v1) const
{
    auto cs = std::make_unique<CoordinateSequence>();
    cs->add(m_vertexList[v0], false);
    cs->add(m_vertexList[v1], false);
    return m_factory->createLineString(std::move(cs));
}


std::vector<std::unique_ptr<LineString>>
SegmentGraph::pairsToGeometry(std::set<std::pair<uint32_t, uint32_t>>& pairs)
{
    // std::cout << "pairsToGeometry pairs.size() = " << pairs.size() << std::endl;

    // Merge all paths we have discovered
    LineMerger lm;
    std::vector<std::unique_ptr<LineString>> lines;
    // Convert vertex pairs into two-point segments
    // Because the set has deduped them we should get
    // a unique collection of segments
    for (const std::pair<uint32_t, uint32_t>& p : pairs) {
        auto ls = pathToGeometry(p.first, p.second);
        lines.emplace_back(ls.release());
    }

    for (const auto& ls : lines) {
        lm.add(ls.get());
    }

    // Stringing together noded collections of unique
    // segments/lines is what LineMerger does best
    return lm.getMergedLineStrings();
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
        
