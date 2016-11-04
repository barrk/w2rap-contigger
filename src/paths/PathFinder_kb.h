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



//using BigDict = HashSet<BigKMer<200>,typename BigKMer<200>::hasher>;

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
        // need to run big kmer pather with lmp reads and pe hbv
        // so actually do just map the reads, the
        hbv.ToLeft(mToLeft);
        hbv.ToRight(mToRight);
        //auto peEdges = hbv.Edges();
        //vecbvec peEdgesbvec; // really should not have to do this like this
        //for (auto edge: peEdges) {peEdgesbvec.push_back(edge)};
        //BigKEdgeBuilder<200>::buildEdges(peDict,&peEdgesbvec);

    }


    //Graph-related methods
    //std::vector<std::vector<uint64_t>> AllPathsFromTo(std::vector<uint64_t> in_edges, std::vector<uint64_t> out_edges, uint64_t max_length);


    //ReadPath-related methods


    //void classify_forks();//how many forks of each type are there?
    //std::vector<uint64_t> best_path_fw(uint64_t edge, int distance); //finds the best path forward for an edge
    //std::array<uint64_t,3> transition_votes(uint64_t left_e,uint64_t right_e);
    //std::array<uint64_t,3> path_votes(std::vector<uint64_t> path);
    //std::array<uint64_t,3> multi_path_votes(std::vector<std::vector<uint64_t>> path);
    //bool path_absolute_best(std::vector<uint64_t> path); //checks a path and its reverse, checks alternatives, true if shold be replaced
    //void untangle_path(std::vector<uint64_t> path);
    //void untangle_pins();
    //void unroll_loops(uint64_t min_side_sizes);//untangles all single choices when support is uncontested
    void untangle_complex_in_out_choices(uint64_t large_frontier_size, bool verbose_separation=false);
    void init_prev_next_vectors();
    //std::vector<std::vector<uint64_t>> is_unrollable_loop(uint64_t e,uint64_t min_side_sizes);//returns size of the unrolled loop
    //uint64_t paths_per_kbp(uint64_t e);
    //std::string edge_pstr(uint64_t e);
    std::string path_str(std::vector<uint64_t> e);
    //std::map<uint64_t,std::vector<uint64_t>> separate_path(std::vector<uint64_t> p, bool verbose_separation=false);
    //bool join_edges_in_path(std::vector<uint64_t> p);
    std::array<std::vector<uint64_t>,2>  get_all_long_frontiers(uint64_t e,uint64_t large_frontier_size);
    //void migrate_readpaths(std::map<uint64_t,std::vector<uint64_t>> edgemap);
    void mapEdgesToLMPReads();


private:
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