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

#include <geos/export.h>

#include <memory>

// Forward declarations
namespace geos {
namespace geom {
class Geometry;
class MultiLineString;
class Polygon;
}
namespace operation {
namespace skeletonize {
}
}
}

namespace geos {
namespace operation { // geos::operation
namespace skeletonize { // geos::operation::skeletonize

class GEOS_DLL Skeletonizer {


public:

    Skeletonizer(const geos::geom::Polygon &poly)
        : inputPolygon(poly)
        {};

    std::unique_ptr<geom::MultiLineString> skeletonize();

    static std::unique_ptr<geom::MultiLineString> skeletonize(const geom::Polygon &poly);

private:

    const geom::Polygon &inputPolygon;

};


} // namespace geos::operation::skeletonize
} // namespace geos::operation
} // namespace geos
