//
// Created by Katie Barr (EI) on 04/11/2016.
//


#include <paths/HyperBasevector.h>
#include <kmers/kmatch/KMatch.h>
#include "lmp_mapper.h"


LMPMapper::LMPMapper(vecbvec lmp_reads, HyperBasevector& hbv, KMatch kmatch): lmp_reads(lmp_reads), hbv(hbv), kMatch(kmatch), read_edge_maps(initalise_read_edge_map()) {}

void LMPMapper::mapReads(){
    kMatch.Hbv2Map(&hbv);
    std::cout << kMatch.edgeMap.size() << std::endl;
    //std::vector<edgeKmerPosition> res = kmatch.lookupRead(lmp_data[0].ToString());
    for (int i=0; i < lmp_reads.size(); i++){
        read_edge_maps.push_back(kMatch.lookupRead(lmp_reads[i].ToString()));
    }
}

// now need a load of helpers- get consecutive offsts from edge vector

void LMPMapper::findFullyMappedEdges(){
    //for (std::vector<edgeKmerPosition>::iterator it = read_edge_maps.begin(); it != read_edge_maps.end(); ++it){
    for (auto read_mapping: read_edge_maps){
        //std::vector<int> getFullyMappedEdges(read_edge_maps[1]);
        std::vector<int> res = getFullyMappedEdges(read_mapping);
}
}


// assume perfect mapping to start with
std::vector<int> LMPMapper::getFullyMappedEdges(std::vector<edgeKmerPosition> read_mapping, int k=31){
    int current_edge_id = -1;
    int consecutive_offsets = 0;
    int last_offset = read_mapping[0].offset;
    std::vector<int> paths;
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
            if (consecutive_offsets == hbv.EdgeObject(it->edge_id).size() - k + 1) {
                paths.push_back(current_edge_id);
            }
            current_edge_id = it->edge_id;
            consecutive_offsets = 0;
        }
        last_offset = it->offset;
        // if the above conditional would miss this because we're at the end of the loop, exit here
        if (std::next(it) == read_mapping.end() && consecutive_offsets == hbv.EdgeObject(it->edge_id).size() - k + 1){
            paths.push_back(current_edge_id);
        }
    }
    return paths;
}
