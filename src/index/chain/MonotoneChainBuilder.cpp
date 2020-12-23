/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: index/chain/MonotoneChainBuilder.java r388 (JTS-1.12)
 *
 **********************************************************************/

#include <geos/index/chain/MonotoneChainBuilder.h>
#include <geos/index/chain/MonotoneChain.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/Quadrant.h>

#include <cassert>
#include <cstdio>
#include <vector>

#ifndef GEOS_DEBUG
#define GEOS_DEBUG 0
#endif

#if GEOS_DEBUG
#include <iostream>
#endif


using namespace geos::geom;

namespace geos {
namespace index { // geos.index
namespace chain { // geos.index.chain


/* static public */
void
MonotoneChainBuilder::getChains(const CoordinateSequence* pts, const void* context,
                                std::vector<MonotoneChain>& mcList)
{
    std::size_t chainStart = 0;
    do {
        std::size_t chainEnd = findChainEnd(*pts, chainStart);
        mcList.emplace_back(*pts, chainStart, chainEnd, context);
        chainStart = chainEnd;
    }
    while (chainStart < (pts->size() - 1));
}


/* private static */
std::size_t
MonotoneChainBuilder::findChainEnd(const CoordinateSequence& pts, std::size_t start)
{

    const std::size_t npts = pts.getSize(); // cache

    assert(start < npts);
    assert(npts); // should be implied by the assertion above,
    // 'start' being unsigned

    std::size_t safeStart = start;

    // skip any zero-length segments at the start of the sequence
    // (since they cannot be used to establish a quadrant)
    while(safeStart < npts - 1
            && pts[safeStart].equals2D(pts[safeStart + 1])) {
        ++safeStart;
    }

    // check if there are NO non-zero-length segments
    if(safeStart >= npts - 1) {
        return npts - 1;
    }

    // determine overall quadrant for chain
    // (which is the starting quadrant)
    int chainQuad = Quadrant::quadrant(pts[safeStart],
                                       pts[safeStart + 1]);

    const Coordinate* prev; // avoid repeated coordinate access by index (virtual call)
    const Coordinate* curr = &pts[start];

    for(std::size_t last = start + 1; last < npts; last++) {
        prev = curr;
        curr = &pts[last];

        // skip zero-length segments, but include them in the chain
        if(!prev->equals2D(*curr)) {
            // compute quadrant for next possible segment in chain
            int quad = Quadrant::quadrant(*prev, *curr);
            if(quad != chainQuad) {
                return last - 1;
            }
        }
    }
#if GEOS_DEBUG
    std::cerr << "MonotoneChainBuilder::findChainEnd() returning" << std::endl;
#endif

    return npts - 1;
}

} // namespace geos.index.chain
} // namespace geos.index
} // namespace geos

