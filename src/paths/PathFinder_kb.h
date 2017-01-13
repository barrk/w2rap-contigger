//
// Created by Katie Barr (EI) on 01/11/2016.
//

#ifndef W2RAP_CONTIGGER_PATHFINDER_H
#define W2RAP_CONTIGGER_PATHFINDER_H

#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "feudal/HashSet.h"
#include "kmers/BigKPather.h"
#include "kmers/kmatch/KMatch.h"
#include "paths/long/large/GapToyTools.h"
#include "kmers/BigKMer.h"
#include "Vec.h"
#include "GFADump.h"
#include <lmp/lmp_mapper.h>
#include <paths/ComplexRegion.h>
#include <fstream>
#include <algorithm>
#include <iostream>


class PathFinderkb {
public:


    PathFinderkb( HyperBasevector& hbv, vec<int>& inv, ReadPathVec& paths, VecULongVec& invPaths,vecbvec& lmp_data, int min_reads = 5) :
            mHBV(hbv),
            mInv(inv),
            mPaths(paths),
            mEdgeToPathIds(invPaths),
            mMinReads(min_reads),
            lmp_data(lmp_data)
    {
        hbv.ToLeft(mToLeft);
        hbv.ToRight(mToRight);
    }


    //ReadPath-related methods
    void init_prev_next_vectors();
    std::string path_str(std::vector<uint64_t> e);
    std::string read_path_str(ReadPath path);
    std::tuple <std::vector<LMPPair >, std::vector<LMPPair >, std::map<uint64_t, std::vector<int> > > mapEdgesToLMPReads();//(std::vector<LMPPair > & lmp_pairs_for_scaffolding);
    void migrate_readpaths(std::map<uint64_t,std::vector<uint64_t>> edgemap);
    std::map<uint64_t,std::vector<uint64_t>> separate_path(std::vector<uint64_t> p, bool verbose_separation);
    void resolveComplexRegionsUsingLMPData();

private:
    void edges_beyond_distance(std::vector<uint64_t>  & long_fronteirs, std::vector<std::vector<uint64_t> >  & paths_to_long_fronteirs, std::vector<uint64_t> &  intermediate_path, uint64_t e, vector<uint64_t > & traversed_edge_list, uint64_t large_frontier_size, int recursion_depth=0, int distance_traversed=0, std::string direction="right");
    HyperBasevector& mHBV;
    vecbvec& lmp_data;
    vec<int>& mInv;
    ReadPathVec& mPaths;
    VecULongVec& mEdgeToPathIds;
    vec<int> mToLeft;
    vec<int> mToRight;
    std::vector<std::vector<uint64_t>> next_edges,prev_edges;
    int mMinReads;


};


#endif //W2RAP_CONTIGGER_PATHFINDER_H