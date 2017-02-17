///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Various things go here, including all abnormal termination error messages for
// Discovar (not including internal errors).  This allows us to maintain some
// consistency between these messages.

#ifndef DISCOVAR_TOOLS_H
#define DISCOVAR_TOOLS_H

#include "CoreTools.h"
#include "Fastavector.h"
#include "paths/HyperBasevector.h"
#include "paths/long/Logging.h"

namespace DiscovarTools{

    // Exit 0 upon finding that the assembly is empty, or there are no reads.

    void ExitAssemblyEmpty( );
    void ExitPathsEmpty( );
    void ExitNoReads( );
    void ExitNoCorrectedReads( );
    void ExitShortReads( const String& additional_info="" );

    // Final message and exit 1 when something really wrong.

    void DiscovarUnhappy( );


    //check if parsed REGIONS input is valid or not
    void CheckDiscovarRegions( const String& REGIONS );

    //check if parsed OUT_HEAD input is valid or not
    void CheckDiscovarOutHead( const String& OUT_HEAD );

    //check if parsed TMP input is valid or not
    void CheckDiscovarTmp( const String& TMP );

    //check if parsed READS input is valid or not
    void CheckDiscovarReads( const String& READS );

    void TestDiscovarRegionsBamsCompatibility( const String& REGIONS,
         const vec<String>& bams );

    void ExitSamtoolsFailed( );

    void CheckReferenceInput( const String& REFERENCE, const String& OUT_HEAD );

    // Forbids region longer than certain length
    const size_t MaxRegionSize=50000000; // Maximum allowable region size
    void CheckRegionSize(size_t size,const String&sPrefix=""); // helper function: check if 'size' is longer than MaxRegionSize

    // Forbids total BAM file size over certain number of bytes, typically called when REGIONS="all", which is then translated to REGIONS==""
    const size_t MaxBAMSize=10ULL*1024*1024*1024; // Maximum allowable BAM size
    void CheckBAMSize(size_t size,const String& REGIONS); // helper function: check if 'size' is larger than MaxBAMSize if REGIONS=="all" or ""

}//namespace DiscovarTools

void SetThreads( uint& NUM_THREADS, const Bool malloc_per_thread_check = True );

#endif
