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

    // Coordinate-to-node-ID mapping (same pattern as SpanningTree)
    struct NodeMapping {
        std::unordered_map<CoordinateXY, std::size_t, CoordinateXY::HashCode> coordToId;

        std::size_t getId(const CoordinateXY& c) {
            auto it = coordToId.find(c);
            if (it != coordToId.end()) return it->second;
            std::size_t id = coordToId.size();
            coordToId[c] = id;
            return id;
        }

        std::size_t size() const {
            return coordToId.size();
        }

        // Returns nodeId → CoordinateXY vector (for snap and centroid)
        std::vector<CoordinateXY> toVector() const {
            std::vector<CoordinateXY> v(coordToId.size());
            for (const auto& kv : coordToId)
                v[kv.second] = kv.first;
            return v;
        }
    };

    // Adjacency list edge
    struct SPEdge {
        std::size_t to;
        double weight;
        std::size_t originalIndex;
    };

    // Result of one Dijkstra run
    struct DijkstraResult {
        std::vector<double>      dist;     // shortest distance from source
        std::vector<std::size_t> pred;     // predecessor node (-1 = none)
        std::vector<std::size_t> predEdge; // curve index used to reach node
    };

    // Priority queue entry for min-heap Dijkstra
    struct PQEntry {
        double dist;
        std::size_t node;
        bool operator>(const PQEntry& o) const { return dist > o.dist; }
    };

    // Build graph from curves (two-pass: register nodes, then add edges)
    void buildGraph(
        const std::vector<const Curve*>& curves,
        NodeMapping& mapping,
        std::vector<std::vector<SPEdge>>& adj)
    {
        // Pass 1: register all nodes
        for (std::size_t i = 0; i < curves.size(); ++i) {
            const Curve* c = curves[i];
            if (!c || c->isEmpty()) continue;
            auto sp = c->getStartPoint();
            auto ep = c->getEndPoint();
            if (!sp || !ep) continue;
            CoordinateXY s(*(sp->getCoordinate()));
            CoordinateXY e(*(ep->getCoordinate()));
            if (s.equals2D(e)) continue;
            mapping.getId(s);
            mapping.getId(e);
        }

        adj.resize(mapping.size());

        // Pass 2: add directed edges in both directions (undirected graph)
        for (std::size_t i = 0; i < curves.size(); ++i) {
            const Curve* c = curves[i];
            if (!c || c->isEmpty()) continue;
            auto sp = c->getStartPoint();
            auto ep = c->getEndPoint();
            if (!sp || !ep) continue;
            CoordinateXY s(*(sp->getCoordinate()));
            CoordinateXY e(*(ep->getCoordinate()));
            if (s.equals2D(e)) continue;
            std::size_t u = mapping.getId(s);
            std::size_t v = mapping.getId(e);
            double len = c->getLength();
            adj[u].push_back({v, len, i});
            adj[v].push_back({u, len, i});
        }
    }

    // Standard lazy Dijkstra from source node
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

            // Stale entry check
            if (d > res.dist[u]) continue;

            for (const SPEdge& e : adj[u]) {
                double nd = res.dist[u] + e.weight;
                if (nd < res.dist[e.to]) {
                    res.dist[e.to] = nd;
                    res.pred[e.to] = u;
                    res.predEdge[e.to] = e.originalIndex;
                    pq.push({nd, e.to});
                }
            }
        }
        return res;
    }

    // Snap a coordinate to the nearest node in the mapping (Euclidean)
    std::size_t snapToNode(
        const CoordinateXY& pt,
        const NodeMapping& mapping)
    {
        const auto coords = mapping.toVector();
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

    // Walk predecessor map from endNode back to startNode, collecting edge indices.
    // Returns edges ordered start → end (position 1 = first edge from start).
    // Returns empty vector if endNode is unreachable.
    std::vector<std::size_t> extractPath(
        std::size_t startNode,
        std::size_t endNode,
        const DijkstraResult& dr)
    {
        if (dr.dist[endNode] == std::numeric_limits<double>::infinity())
            return {};

        std::vector<std::size_t> edgesReverse;
        std::size_t cur = endNode;
        while (cur != startNode) {
            std::size_t edgeIdx = dr.predEdge[cur];
            if (edgeIdx == std::size_t(-1)) return {}; // disconnected
            edgesReverse.push_back(edgeIdx);
            cur = dr.pred[cur];
        }
        std::reverse(edgesReverse.begin(), edgesReverse.end());
        return edgesReverse;
    }

    // Return the index of the node with maximum finite distance, or
    // std::size_t(-1) if all distances are infinite.
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
    result.assign(curves.size(), 0);
    if (curves.empty()) return;

    NodeMapping mapping;
    std::vector<std::vector<SPEdge>> adj;
    buildGraph(curves, mapping, adj);

    if (mapping.size() == 0) return;

    std::size_t startNode = snapToNode(start, mapping);
    std::size_t endNode   = snapToNode(end,   mapping);

    if (startNode == endNode) return;

    DijkstraResult dr = runDijkstra(startNode, adj);

    std::vector<std::size_t> path = extractPath(startNode, endNode, dr);

    for (std::size_t pos = 0; pos < path.size(); ++pos) {
        result[path[pos]] = pos + 1;
    }
}

void
ShortestPath::longestShortestPath(
    const std::vector<const geom::Curve*>& curves,
    std::vector<std::size_t>& result)
{
    result.assign(curves.size(), 0);
    if (curves.empty()) return;

    NodeMapping mapping;
    std::vector<std::vector<SPEdge>> adj;
    buildGraph(curves, mapping, adj);

    if (mapping.size() == 0) return;

    const std::vector<CoordinateXY> coords = mapping.toVector();
    const std::size_t n = coords.size();

    // Step 1: Find centroid of all node coordinates
    double cx = 0.0, cy = 0.0;
    for (const auto& c : coords) { cx += c.x; cy += c.y; }
    cx /= static_cast<double>(n);
    cy /= static_cast<double>(n);
    CoordinateXY centroid(cx, cy);

    // Step 2: Snap centroid to nearest node → "center"
    std::size_t center = snapToNode(centroid, mapping);

    // Step 3: Dijkstra from center → find node A (furthest by graph distance)
    DijkstraResult drCenter = runDijkstra(center, adj);
    std::size_t nodeA = furthestReachableNode(drCenter.dist);
    if (nodeA == std::size_t(-1)) return; // no reachable nodes

    // Step 4: Dijkstra from A → find node B and save predecessor map for A→C
    DijkstraResult drA = runDijkstra(nodeA, adj);
    std::size_t nodeB = furthestReachableNode(drA.dist);
    if (nodeB == std::size_t(-1)) return;

    // Step 5: Dijkstra from B → find node C (furthest from B)
    DijkstraResult drB = runDijkstra(nodeB, adj);
    std::size_t nodeC = furthestReachableNode(drB.dist);
    if (nodeC == std::size_t(-1)) return;

    // Step 6: Extract path A→C using predecessor map from drA.
    // If nodeC == nodeA the graph is symmetric (chain, star, etc.) and the
    // triple-sweep degenerates; fall back to the A→B diameter path.
    std::size_t target = (nodeC != nodeA) ? nodeC : nodeB;
    std::vector<std::size_t> path = extractPath(nodeA, target, drA);

    // Step 7: Encode result with position (1-indexed)
    for (std::size_t pos = 0; pos < path.size(); ++pos) {
        result[path[pos]] = pos + 1;
    }
}

} // namespace geos::operation::shortest
} // namespace geos::operation
} // namespace geos
