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

/*
 * when a region which is solveable is found, create an object for it to hold the paths
 * then, consistency between regions can be validated before tying to separate,
 * and paths which resolve a regin can be selected over ones that don't
 */

namespace LMPRegions {
    class ComplexRegion {
    public:
        ComplexRegion::ComplexRegion(std::vector<uint64_t> edges_in, std::vector<uint64_t> edges_out,
                                     vec<int> &involution, int insert_size = 5000);

        ComplexRegion::ComplexRegion();

        std::vector<uint64_t> edges_in;
        std::vector<uint64_t> edges_out;
        std::vector<std::vector<uint64_t> > paths;
        vec<int> involution;

        void AddPath(ReadPath path);

        bool solved;
        int insert_size;
    private:
        std::vector<uint64_t> caonicalisePath(ReadPath &path);

        std::vector<std::vector<uint64_t> > candidate_paths;
        std::vector<std::vector<uint64_t> > selected_paths;
        // to ensure all paths are on same part of the graph, in same direction
        std::vector<uint64_t> edges_in_canonical;
        std::vector<uint64_t> edges_out_canonical;

        std::vector<uint64_t> ComplexRegion::canonicalisePath(ReadPath path);


    };

/*
class ComplexRegionCollection {
public:
    ComplexRegionCollection::ComplexRegionCollection();
    void AddRegion(ComplexRegion complex_region);
    bool ContainsRegionWithEdgesIn(std::vector<uint64_t>);
    ComplexRegion GetRegionWithEdgesIn(std::vector<uint64_t>);
    std::vector<ComplexRegion> complex_regions;

private:
    std::map<std::vector<uint64_t>, int> edges_in_to_region_index;
    //bool OverlapsOtherRegions();


};
*/
}

#endif //W2RAP_CONTIGGER_COMPLEXREGION_H
