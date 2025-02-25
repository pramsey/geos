
# Skeletonizer Plan

## Pre-Processing

- **DONE - Densification** : In the case of a ring with start/end point at an acute angle, the densification can lead to bad results (edges outside polygon). Need to detect case and move start/end to safe, oblique angle.

- **Bad Inputs** : Remove any zero-length and small-length segments

- **Bad Inputs** : DP away unnecessary vertices before adding in the new dense points. Hopefully removes very short edges. Have to do so without losing the input/output vertices.

- **Bad Inputs** : Remove spikes, gores, and zigzags?

- **Bad Inputs** : Search for edges that are spatially "close" while also being a long distance "along the ring" from the query edge. These are places where the ring folds in on itself, or has a narrow channel with another ring. Densify these edges even more. Build STRee, and test every edge for this condition.

- **Densification** : Calculate correct level of densification? Relative to object size? Relative to smallest segment? Relative to minimum clearance?

- **Input/Output** : Replace input/output vertices with two very close together points to drive the diagram through the correct position in the point field.

## Post-Voronoi

- How to remove unnecessary edges? Full containment as the condition might remove edges required to connect input/output points to the central spine.

- Build collection of edges output into a graph and perform the "furthest point" calculation to figure out longest interior skeleton line

- Deal with case in which there are multiple disconnected graphs. Find the largest, ignore the small ones?
