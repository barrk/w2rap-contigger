///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLE_GAPS_H
#define ASSEMBLE_GAPS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"

void AssembleGaps2( HyperBasevector& hb, vec<int>& inv2, ReadPathVec& paths2, 
     const vecbasevector& bases, VecPQVec const& quals,
     const String& work_dir, std::vector<int>,
     vecbvec& new_stuff, const Bool CYCLIC_SAVE,
     const int A2V, const int MAX_PROX_LEFT,
     const int MAX_PROX_RIGHT, const int max_copies, const int max_bpaths, const int max_iterations,
     const int pair_sample, const std::string out_lr_dir, const std::string in_lr_dir, const bool dump_lr );

#endif
