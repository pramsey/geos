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

#include <geos/operation/shortest/ShortestPath.h>
#include <geos/geom/Curve.h>
#include <geos/geom/Point.h>
#include <geos/geom/Coordinate.h>

#include <vector>
#include <unordered_map>
#include <queue>
#include <limits>
#include <algorithm>
#include <numeric>

using namespace geos::geom;

namespace geos {
namespace operation { // geos::operation
namespace shortest { // geos::operation::shortest

namespace {

    /**
     * \brief Coordinate-to-node-ID mapping.
     *
     * Maps unique CoordinateXY values to sequential size_t IDs (0 to N-1).
     * Maintains an inverse mapping (idToCoord) for efficient coordinate lookup.
     */
    struct NodeMapping {
        std::unordered_map<CoordinateXY, std::size_t, CoordinateXY::HashCode> coordToId;
        std::vector<CoordinateXY> idToCoord;

        /** \brief Returns the ID for a coordinate, creating it if necessary. */
        std::size_t getId(const CoordinateXY& c) {
            auto it = coordToId.find(c);
            if (it != coordToId.end()) return it->second;
            std::size_t id = idToCoord.size();
            coordToId[c] = id;
            idToCoord.push_back(c);
            return id;
        }

        /** \brief Returns the number of unique nodes registered. */
        std::size_t size() const {
            return idToCoord.size();
        }

        /** \brief Returns the mapping from ID to CoordinateXY. */
        const std::vector<CoordinateXY>& toVector() const {
            return idToCoord;
        }
    };

    /** \brief Adjacency list edge representation. */
    struct SPEdge {
        std::size_t to;             // Target node ID
        double weight;              // Edge weight (curve length)
        std::size_t originalIndex;  // Index in the input curves vector
    };

    /** \brief Result of a Dijkstra algorithm execution. */
    struct DijkstraResult {
        std::vector<double>      dist;     // Shortest distance from source to each node
        std::vector<std::size_t> pred;     // Predecessor node ID (-1 if unreachable/none)
        std::vector<std::size_t> predEdge; // Original curve index used to reach the node
    };

    /** \brief Priority queue entry for min-heap Dijkstra. */
    struct PQEntry {
        double dist;
        std::size_t node;
        // Priority queue is a max-heap by default; operator> makes it a min-heap.
        bool operator>(const PQEntry& o) const { return dist > o.dist; }
    };

    /**
     * \brief Build an undirected graph from a vector of Curves.
     *
     * Nodes are formed from the unique start and end points of the curves.
     * Edges are added in both directions. Closed curves (rings) are currently skipped.
     */
    void buildGraph(
        const std::vector<const Curve*>& curves,
        NodeMapping& mapping,
        std::vector<std::vector<SPEdge>>& adj)
    {
        // Temporary storage for edges before we know the final node count.
        struct TempEdge { std::size_t u, v; double len; std::size_t idx; };
        std::vector<TempEdge> tempEdges;
        tempEdges.reserve(curves.size());

        for (std::size_t i = 0; i < curves.size(); ++i) {
            const Curve* c = curves[i];
            if (!c || c->isEmpty()) continue;
            auto sp = c->getStartPoint();
            auto ep = c->getEndPoint();
            if (!sp || !ep) continue;
            CoordinateXY s(*(sp->getCoordinate()));
            CoordinateXY e(*(ep->getCoordinate()));
            
            // Skip self-loops (curves starting and ending at the same node)
            if (s.equals2D(e)) continue;

            std::size_t u = mapping.getId(s);
            std::size_t v = mapping.getId(e);
            tempEdges.push_back({u, v, c->getLength(), i});
        }

        // Initialize adjacency list with the correct number of nodes
        adj.assign(mapping.size(), {});
        for (const auto& te : tempEdges) {
            // Undirected graph: add edge in both directions
            adj[te.u].push_back({te.v, te.len, te.idx});
            adj[te.v].push_back({te.u, te.len, te.idx});
        }
    }

    /**
     * \brief Standard Dijkstra's algorithm using a priority queue (min-heap).
     *
     * Computes the shortest path distance and predecessors from a single source node.
     */
    DijkstraResult runDijkstra(
        std::size_t source,
        const std::vector<std::vector<SPEdge>>& adj)
    {
        const std::size_t n = adj.size();
        DijkstraResult res;
        res.dist.assign(n, std::numeric_limits<double>::infinity());
        res.pred.assign(n, std::size_t(-1));
        res.predEdge.assign(n, std::size_t(-1));

        res.dist[source] = 0.0;

        using PQ = std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>>;
        PQ pq;
        pq.push({0.0, source});

        while (!pq.empty()) {
            PQEntry top = pq.top();
            pq.pop();
            double d = top.dist;
            std::size_t u = top.node;

            // If we found a shorter path already, this is a stale entry.
            if (d > res.dist[u]) continue;

            for (const SPEdge& e : adj[u]) {
                double newDist = res.dist[u] + e.weight;
                if (newDist < res.dist[e.to]) {
                    res.dist[e.to] = newDist;
                    res.pred[e.to] = u;
                    res.predEdge[e.to] = e.originalIndex;
                    pq.push({newDist, e.to});
                }
            }
        }
        return res;
    }

    /** \brief Snap a coordinate to the nearest graph node using squared Euclidean distance. */
    std::size_t snapToNode(
        const CoordinateXY& pt,
        const std::vector<CoordinateXY>& coords)
    {
        if (coords.empty()) return std::size_t(-1);

        std::size_t bestId = 0;
        double bestDist = pt.distanceSquared(coords[0]);
        for (std::size_t i = 1; i < coords.size(); ++i) {
            double d = pt.distanceSquared(coords[i]);
            if (d < bestDist) {
                bestDist = d;
                bestId = i;
            }
        }
        return bestId;
    }

    /**
     * \brief Trace back from endNode to startNode using the predecessor map.
     *
     * Returns a vector of curve indices forming the path, in order from start to end.
     */
    std::vector<std::size_t> extractPath(
        std::size_t startNode,
        std::size_t endNode,
        const DijkstraResult& dr)
    {
        // If endNode was never reached, return an empty path.
        if (dr.dist[endNode] == std::numeric_limits<double>::infinity())
            return {};

        std::vector<std::size_t> edgesReverse;
        std::size_t cur = endNode;
        while (cur != startNode) {
            std::size_t edgeIdx = dr.predEdge[cur];
            // If we hit a node with no predecessor before reaching startNode, it's a bug or disconnection.
            if (edgeIdx == std::size_t(-1)) return {}; 
            edgesReverse.push_back(edgeIdx);
            cur = dr.pred[cur];
        }
        // Reverse to get path from start to end.
        std::reverse(edgesReverse.begin(), edgesReverse.end());
        return edgesReverse;
    }

    /** \brief Returns the index of the reachable node furthest from the source. */
    std::size_t furthestReachableNode(const std::vector<double>& dist)
    {
        std::size_t best = std::size_t(-1);
        double bestDist = -1.0;
        for (std::size_t i = 0; i < dist.size(); ++i) {
            if (dist[i] != std::numeric_limits<double>::infinity() && dist[i] > bestDist) {
                bestDist = dist[i];
                best = i;
            }
        }
        return best;
    }

} // anonymous namespace

void
ShortestPath::shortestPath(
    const std::vector<const geom::Curve*>& curves,
    const geom::CoordinateXY& start,
    const geom::CoordinateXY& end,
    std::vector<std::size_t>& result)
{
    // Initialize result with zeros (size equal to number of input curves)
    result.assign(curves.size(), 0);
    if (curves.empty()) return;

    // 1. Build the graph adjacency list from input curves
    NodeMapping mapping;
    std::vector<std::vector<SPEdge>> adj;
    buildGraph(curves, mapping, adj);

    if (mapping.size() == 0) return;

    // 2. Snap input start/end points to the nearest graph nodes
    const auto& coords = mapping.toVector();
    std::size_t startNode = snapToNode(start, coords);
    std::size_t endNode   = snapToNode(end,   coords);

    // If start and end snap to the same node, the shortest path is empty.
    if (startNode == endNode) return;

    // 3. Run Dijkstra's algorithm from the start node
    DijkstraResult dr = runDijkstra(startNode, adj);

    // 4. Extract the sequence of curve indices forming the path
    std::vector<std::size_t> path = extractPath(startNode, endNode, dr);

    // 5. Populate the result vector with 1-based positions
    for (std::size_t pos = 0; pos < path.size(); ++pos) {
        result[path[pos]] = pos + 1;
    }
}

void
ShortestPath::longestShortestPath(
    const std::vector<const geom::Curve*>& curves,
    std::vector<std::size_t>& result)
{
    // Initialize result with zeros
    result.assign(curves.size(), 0);
    if (curves.empty()) return;

    // 1. Build the graph
    NodeMapping mapping;
    std::vector<std::vector<SPEdge>> adj;
    buildGraph(curves, mapping, adj);

    if (mapping.size() == 0) return;

    const std::size_t nNodes = mapping.size();
    const auto& allCoords = mapping.toVector();
    std::vector<bool> visited(nNodes, false);
    
    double maxDiameter = -1.0;
    std::vector<std::size_t> bestPath;

    // 2. Iterate through all connected components to find the global diameter.
    for (std::size_t i = 0; i < nNodes; ++i) {
        if (visited[i]) continue;

        // A. Identify all nodes in the current connected component using BFS.
        std::vector<std::size_t> componentNodes;
        std::queue<std::size_t> q;
        q.push(i);
        visited[i] = true;
        while (!q.empty()) {
            std::size_t u = q.front(); q.pop();
            componentNodes.push_back(u);
            for (const auto& e : adj[u]) {
                if (!visited[e.to]) {
                    visited[e.to] = true;
                    q.push(e.to);
                }
            }
        }

        // B. Heuristic for finding component diameter: Triple-Dijkstra Sweep.
        
        // i. Find the centroid of nodes in this component and snap it to a "center" node.
        double cx = 0.0, cy = 0.0;
        for (std::size_t nodeIdx : componentNodes) {
            cx += allCoords[nodeIdx].x;
            cy += allCoords[nodeIdx].y;
        }
        cx /= static_cast<double>(componentNodes.size());
        cy /= static_cast<double>(componentNodes.size());
        CoordinateXY centroid(cx, cy);

        std::size_t center = componentNodes[0];
        double bDist = centroid.distanceSquared(allCoords[center]);
        for (std::size_t j = 1; j < componentNodes.size(); ++j) {
            double d = centroid.distanceSquared(allCoords[componentNodes[j]]);
            if (d < bDist) {
                bDist = d;
                center = componentNodes[j];
            }
        }

        // ii. First Dijkstra: find node A furthest from the "center".
        DijkstraResult drCenter = runDijkstra(center, adj);
        std::size_t nodeA = furthestReachableNode(drCenter.dist);
        if (nodeA == std::size_t(-1)) continue;

        // iii. Second Dijkstra: find node B furthest from A.
        DijkstraResult drA = runDijkstra(nodeA, adj);
        std::size_t nodeB = furthestReachableNode(drA.dist);
        if (nodeB == std::size_t(-1)) continue;

        // iv. Third Dijkstra: find node C furthest from B.
        DijkstraResult drB = runDijkstra(nodeB, adj);
        std::size_t nodeC = furthestReachableNode(drB.dist);
        if (nodeC == std::size_t(-1)) continue;

        // C. Evaluate the diameter found in this component.
        double diamA = drA.dist[nodeB]; // distance(A, B)
        double diamB = drB.dist[nodeC]; // distance(B, C)
        
        if (std::max(diamA, diamB) > maxDiameter) {
            maxDiameter = std::max(diamA, diamB);
            // If B-C is longer and C is not A, it's a better candidate for diameter.
            if (diamB > diamA && nodeC != nodeA) {
                bestPath = extractPath(nodeB, nodeC, drB);
            } else {
                bestPath = extractPath(nodeA, nodeB, drA);
            }
        }
    }

    // 3. Populate the result vector with the path from the best component found.
    for (std::size_t pos = 0; pos < bestPath.size(); ++pos) {
        result[bestPath[pos]] = pos + 1;
    }
}

} // namespace geos::operation::shortest
} // namespace geos::operation
} // namespace geos
