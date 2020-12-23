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
 * Last port: noding/MCIndexSegmentSetMutualIntersector.java r388 (JTS-1.12)
 *
 **********************************************************************/

#include <geos/noding/MCIndexSegmentSetMutualIntersector.h>
#include <geos/noding/SegmentSetMutualIntersector.h>
#include <geos/noding/SegmentString.h>
#include <geos/noding/SegmentIntersector.h>
#include <geos/index/SpatialIndex.h>
#include <geos/index/chain/MonotoneChain.h>
#include <geos/index/chain/MonotoneChainBuilder.h>
#include <geos/index/chain/MonotoneChainOverlapAction.h>
#include <geos/index/strtree/SimpleSTRtree.h>
// std
#include <cstddef>

using namespace geos::index::chain;

namespace geos {
namespace noding { // geos::noding


/*private*/
void
MCIndexSegmentSetMutualIntersector::addToIndex(const SegmentString* segStr)
{
    MonotoneChainBuilder::getChains(segStr->getCoordinates(),
                                    segStr, indexChains);

}


/*private*/
void
MCIndexSegmentSetMutualIntersector::addToMonoChains(const SegmentString* segStr)
{
    MonotoneChainBuilder::getChains(segStr->getCoordinates(),
                                    segStr, monoChains);

}


/*private*/
void
MCIndexSegmentSetMutualIntersector::intersectChains()
{
    MCIndexSegmentSetMutualIntersector::SegmentOverlapAction overlapAction(*segInt);

    std::vector<void*> overlapChains;
    for(auto& queryChain : monoChains) {
        overlapChains.clear();

        index->query(&(queryChain.getEnvelope()), overlapChains);

        for(std::size_t j = 0, nj = overlapChains.size(); j < nj; j++) {
            MonotoneChain* testChain = (MonotoneChain*)(overlapChains[j]);

            queryChain.computeOverlaps(testChain, &overlapAction);
            nOverlaps++;
            if(segInt->isDone()) {
                return;
            }
        }
    }
}


/* public */
void
MCIndexSegmentSetMutualIntersector::setBaseSegments(SegmentString::ConstVect* segStrings)
{
    // NOTE - mloskot: const qualifier is removed silently, dirty.

    for(const SegmentString* css: *segStrings) {
        SegmentString* ss = const_cast<SegmentString*>(css);
        addToIndex(ss);
    }
}

/*public*/
void
MCIndexSegmentSetMutualIntersector::process(SegmentString::ConstVect* segStrings)
{
    if (!indexBuilt) {
        for (auto& mc: indexChains) {
            index->insert(&(mc.getEnvelope()), &mc);
        }
        indexBuilt = true;
    }

    // Reset counters for new inputs
    monoChains.clear();
    processCounter = indexCounter + 1;
    nOverlaps = 0;

    for(const SegmentString* css: *segStrings) {
        SegmentString* ss = const_cast<SegmentString*>(css);
        addToMonoChains(ss);
    }
    intersectChains();
}


/* public */
void
MCIndexSegmentSetMutualIntersector::SegmentOverlapAction::overlap(
    MonotoneChain& mc1, std::size_t start1, MonotoneChain& mc2, std::size_t start2)
{
    SegmentString* ss1 = (SegmentString*)(mc1.getContext());
    SegmentString* ss2 = (SegmentString*)(mc2.getContext());

    si.processIntersections(ss1, start1, ss2, start2);
}

} // geos::noding
} // geos

