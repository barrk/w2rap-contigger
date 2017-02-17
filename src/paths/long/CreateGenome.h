///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef CREATE_GENOME_H
#define CREATE_GENOME_H

#include "Basevector.h"
#include "CoreTools.h"
#include "IntPairVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/Logging.h"
#include "paths/long/DiscovarTools.h"

class ref_data {

     public:

     ref_data( ) : LG(12) { }

     vecbasevector G, G3, G3plus;
     vec<int> Gplus_ext;
     vec<HyperBasevector> GH;
     vec<double> ploidy;
     vec<bool> is_circular;
     VecIntPairVec Glocs, G3locs, G3pluslocs;
     int LG;
};

class RefTraceControl;

void ParseRegions( const String& X, vec<String>& regions );

// note that this function is duplicated and modded to
// CreateGenome_Discovar, as a quick fix of memory/speed-related issue.
// Any changes here should also be made there

void CreateGenome( const String& IN_GENOME, const String& SAMPLE, const String& X,
     const String& HUMAN_CONTROLS, String& X_actual, ref_data& ref,
     const double GENOME_SUB_PERCENT,RefTraceControl& RTCtrl );


void CreateGlocs( const vecbasevector& G, unsigned const LG, VecIntPairVec& Glocs );





#endif
