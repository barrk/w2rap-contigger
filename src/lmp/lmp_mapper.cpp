//
// Created by Katie Barr (EI) on 04/11/2016.
//


#include <paths/HyperBasevector.h>
#include <kmers/kmatch/KMatch.h>
#include "lmp_mapper.h"
#include <algorithm>

LMPMapper::LMPMapper(vecbvec lmp_reads, HyperBasevector& hbv, KMatch kmatch): lmp_reads(lmp_reads), hbv(hbv), kMatch(kmatch), read_edge_maps(initalise_read_edge_map()) {}

bool compareEdgeKmerPositions(const edgeKmerPosition &ekp1, const edgeKmerPosition &ekp2){
    // sort by edge id first, then offset. so a higher edge id with a lower offset would be the greater one
    if (ekp1.edge_id < ekp2.edge_id){
        return true;
    } else if (ekp1.edge_id == ekp2.edge_id && ekp1.offset < ekp2.offset){
        // we're on the same edge, so the one with the lowest offset is lower
        return true;
    } else {
        return false;
    }

}


/*ReadPath LMPMapper::getReadMathPair(int edge_id);{

}

ReadPath LMPMapper::getReadMathPair(int read_index){
    ReadPath pair = read_index%2 == 0 ? read_paths[read_index + 1] : read_paths[read_index -1 1]
    return pair;
}*/


void LMPMapper::mapReads(){
    kMatch.Hbv2Map(&hbv);
    std::cout << kMatch.edgeMap.size() << std::endl;
    //std::vector<edgeKmerPosition> res = kmatch.lookupRead(lmp_data[0].ToString());
    for (int i=0; i < lmp_reads.size(); i++){
        read_edge_maps.push_back(kMatch.lookupRead(lmp_reads[i].ToString()));
    }
}


void LMPMapper::convertLMPPairsToReadPaths(){
    LMPPair lmp_pair;
    //for (std::vector<edgeKmerPosition>::iterator it = read_edge_maps.begin(); it != read_edge_maps.end(); ++it){
    std::vector<LMPPair > read_paths;
    for (int i=0; i < read_edge_maps.size(); ++i){
        std::vector<edgeKmerPosition> read_mapping_p1 = read_edge_maps[i];
        // ensure that full edges will be together, with offsets in increasing order, in read_mapping vector
        lmp_pair.p1 = sortMappingsFindullyMappedEdges(read_mapping_p1);
        //++read_mapping; got 'no matching function call' error for these when i tried to iterate using type inference
        //std::next(read_mapping);
        std::vector<edgeKmerPosition> read_mapping_p2 = read_edge_maps[i+1];
        lmp_pair.p2 = sortMappingsFindullyMappedEdges(read_mapping_p2);
        read_paths.push_back(lmp_pair); // check if no fully mapped edges are found, something is still added to the vector, as we need original indices to find pairs
    }
}

ReadPath LMPMapper::sortMappingsFindullyMappedEdges(std::vector<edgeKmerPosition>  read_mapping){
    std::sort(read_mapping.begin(), read_mapping.end(), compareEdgeKmerPositions);
    return getFullyMappedEdges(read_mapping);
}


// assume perfect mapping to start with
ReadPath LMPMapper::getFullyMappedEdges(std::vector<edgeKmerPosition> read_mapping, int k=31){
    int current_edge_id = read_mapping[0].edge_id;
    int consecutive_offsets = 0;
    int last_offset = read_mapping[0].offset;
    ReadPath paths;
    for (auto it = read_mapping.begin(); it != read_mapping.end(); ++it) {
        // if we're still on the same edge
        if (it->edge_id == current_edge_id) {
            // if the edge offset is consecutive
            if (it->offset == (last_offset + 1)) {
                consecutive_offsets += 1;
            }
        }
            // if we're onto a new edge, then check if we completed the previous edge -
        else {
            // if the number of consecutive offsets is equal to the number of kmers on the read, -1 because we count transitions its fully mapped, L - k
            if (consecutive_offsets == hbv.EdgeObject(current_edge_id).size() - k) {
                paths.push_back(current_edge_id);
            }
            current_edge_id = it->edge_id;
            consecutive_offsets = 0;
        }
        last_offset = it->offset;
        // if the above conditional would miss this because we're at the end of the loop, exit here
        if (std::next(it) == read_mapping.end() && consecutive_offsets == hbv.EdgeObject(current_edge_id).size() - k){
            paths.push_back(current_edge_id);
        }
    }
    return paths;
}


