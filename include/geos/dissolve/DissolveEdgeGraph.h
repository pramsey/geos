/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (c) 2025 Martin Davis
 * Copyright (C) 2025 Paul Ramsey <pramsey@cleverelephant.ca>
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#pragma once

#include <geos/edgegraph/EdgeGraph.h>
#include <geos/export.h>


// Forward declarations
namespace geos {
namespace geom {
    class CoordinateXYZM;
}
namespace edgegraph {
    class HalfEdge;
}
}


namespace geos {      // geos.
namespace dissolve {  // geos.dissolve


class GEOS_DLL DissolveEdgeGraph : public edgegraph::EdgeGraph {

    using CoordinateXYZM = geos::geom::CoordinateXYZM;
    using HalfEdge = geos::edgegraph::HalfEdge;


private:

    std::vector<std::unique_ptr<HalfEdge>> edgeStore;

public:

    HalfEdge* createEdge(const CoordinateXYZM& p0);

};

} // namespace geos.dissolve
} // namespace geos

