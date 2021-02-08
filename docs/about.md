---
layout: page
title: About
nav_order: 2
permalink: /about/
---

# About

GEOS is a C/C++ library for spatial computational geometry of the sort generally used by "geographic information systems" software. GEOS is a core dependency of [PostGIS](http://postgis.net), [QGIS](http://qgis.org), and [Shapely](https://shapely.readthedocs.io/en/stable/project.html).

## Capabilities

Spatial Model and Functions

* **Geometries**: Point, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon, GeometryCollection
* **Predicates**: Intersects, Touches, Disjoint, Crosses, Within, Contains, Overlaps, Equals, Covers
* **Operations**: Union, Distance, Intersection, Symmetric Difference, Convex Hull, Envelope, Buffer, Simplify, Polygon Assembly, Valid, Area, Length,
* Prepared geometries (pre-spatially indexed)
* STR spatial index
* OGC Well Known Text (WKT) and Well Known Binary (WKB) encoders and decoders.

## API Features

* C++ API (will likely change across versions)
* C API (provides long-term ABI stability)
* Thread safety (using the reentrant API)

## License

GEOS is [open source software](https://opensource.com/resources/what-open-source) available under the terms of ​[GNU Lesser General Public License](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html) (LGPL).

