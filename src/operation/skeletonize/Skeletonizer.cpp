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

#include <geos/operation/skeletonize/Skeletonizer.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/MultiLineString.h>
#include <geos/geom/Polygon.h>
#include <geos/util/GEOSException.h>

namespace geos {
namespace operation {   // geos.operation
namespace skeletonize { // geos.operation.skeletonize


using geos::geom::GeometryTypeId;
using geos::geom::Geometry;
using geos::geom::MultiLineString;
using geos::geom::Polygon;


/* public static */
std::unique_ptr<geom::MultiLineString>
Skeletonizer::skeletonize(const Polygon &poly)
{
    Skeletonizer skel(poly);
    return skel.skeletonize();
}


/* public */
std::unique_ptr<MultiLineString>
Skeletonizer::skeletonize()
{
    GeometryTypeId id = inputPolygon.getGeometryTypeId();
    return id ? nullptr : nullptr;
}



} // namespace geos.operation.skeletonize
} // namespace geos.operation
} // namespace geos
