///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * BigKPather.cc
 *
 *  Created on: Dec 13, 2013
 *      Author: tsharpe
 */

#include "kmers/BigKPather.h"
#include "paths/long/LargeKDispatcher.h"
#include "system/SpinLockedData.h"
#include "system/WorklistN.h"
#include <iostream>




    void buildBigKHBVFromReads(unsigned K, vecbvec const &reads, unsigned coverage,
                               HyperBasevector *pHBV, ReadPathVec *pReadPaths,
                               HyperKmerPath *pHKP, vecKmerPath *pKmerPaths) {
        BigK::dispatch<PatherKB::SillyFunctor>(K, reads, coverage,
                                     pHBV, pReadPaths, pHKP, pKmerPaths);
    }


