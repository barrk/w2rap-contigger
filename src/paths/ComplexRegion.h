//
// Created by Katie Barr (EI) on 15/12/2016.
//

#ifndef W2RAP_CONTIGGER_COMPLEXREGION_H
#define W2RAP_CONTIGGER_COMPLEXREGION_H

#include <cstdint>
#include <utility>
#include <vector>
#include <map>
#include <Vec.h>
#include "paths/long/ReadPath.h"

typedef struct {
    uint64_t edge_id;
    uint64_t translated_edge_id;
    std::vector<int> pair_ids;
    std::vector<uint64_t> path_from_center;
    bool into; // true if this is an edge going into the region
} BoundingEdge;

/*
 * when a region which is solveable is found, create an object for it to hold the paths
 * then, consistency between regions can be validated before tying to separate,
 * and paths which resolve a regin can be selected over ones that don't
 */

    class ComplexRegion {
    public:
        ComplexRegion::ComplexRegion(std::vector<uint64_t> edges_in, std::vector<uint64_t> edges_out,
                                     vec<int> &involution, int insert_size = 5000);

        ComplexRegion::ComplexRegion();

        std::vector<uint64_t> edges_in;
        std::vector<uint64_t> edges_out;
        std::vector<std::vector<uint64_t> > paths;
        vec<int> involution;// would be nicer to have the region 'know' its in a collection, and use involution from there

        void AddPath(ReadPath path);

        bool solved;
        int insert_size;
        // to ensure all paths are on same part of the graph, in same direction
        std::vector<uint64_t> edges_in_canonical;
        std::vector<uint64_t> edges_out_canonical;
        void canonicaliseEdgesInOut();
        void AddPairId(int edge_index, int pair_id, bool in);
        void AddPathTo(int edge_index, bool in, std::vector<uint64_t> path_from_center);
        void FindSolveablePairs();
        void isSolved(int min_count);

    private:
        std::vector<uint64_t> caonicalisePath(ReadPath &path);

        std::vector<std::vector<uint64_t> > candidate_paths;
        std::vector<std::vector<uint64_t> > selected_paths;

        std::vector<uint64_t> canonicalisePath(ReadPath path);

        std::map<uint64_t, std::vector<int> > pair_ids;
        std::map<uint64_t, uint64_t>  edge_translations;
        std::vector<BoundingEdge> edges_in_detailed;
        std::vector<BoundingEdge> edges_out_detailed;
        std::map<std::pair<uint64_t, uint64_t>, int > combination_counts;

        //void ComplexRegion::canonicaliseEdgesInOut(std::vector<uint64_t> edges_in, std::vector<uint64_t> edges_out);

    };


class ComplexRegionCollection {
public:
    ComplexRegionCollection(vec<int> &involution);
    vec<int> involution;
    //void AddRegion(ComplexRegion complex_region);
    void AddRegion(std::vector<uint64_t> edges_in, std::vector<uint64_t> edges_out,
                   vec<int> &involution, int insert_size = 5000);
    bool ContainsRegionWithEdges(std::vector<uint64_t> edges_in, std::vector<uint64_t> edges_out);
    ComplexRegion GetRegionWithEdges(std::vector<uint64_t> edges_in, std::vector<uint64_t> edges_out);
    std::vector<ComplexRegion> complex_regions;

private:
    std::map<std::pair< std::vector<uint64_t>, std::vector<uint64_t> >, int> edges_to_region_index;
    std::pair< std::vector<uint64_t>, std::vector<uint64_t> > canonicaliseEdgesInOut(std::vector<uint64_t> edges_in, std::vector<uint64_t> edges_out);
    //bool OverlapsOtherRegions();


};



#endif //W2RAP_CONTIGGER_COMPLEXREGION_H
