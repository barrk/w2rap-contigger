//
// Created by Katie Barr (EI) and Rob Vickerstaff (EMR) on 24/02/2017.
//

#include <kmers/kmatch/KMatch.h>
#include "marker_mapping.h"

MarkerMapper::MarkerMapper(vecbvec& marker_sequences, std::vector<std::string> marker_ids, HyperBasevector& hbv, vec<int>& inv, KMatch kmatch): marker_sequences(marker_sequences), marker_ids(marker_ids), hbv(hbv), inv(inv) ,kMatch(kmatch) {
    hbv.ToLeft(mToLeft);
    hbv.ToRight(mToRight);
}



void MarkerMapper::mapMarkers(){
    kMatch.Hbv2Map(&hbv);
    //PrintMemUsage(); // 1125576
    std::cout << Date() << "starting marker probe sequence mapping"<< std::endl;
    int mapped_to_single_edge = 0;
    int mapped_to_multiple_edge = 0;
    int unmapped_marker = 0;
    std::map<uint64_t, std::atomic<int>> mapping_counts;
    for (int i=0; i < marker_sequences.size(); i++){
        std::string read = (marker_sequences)[i].ToString();
        // make guess that if mapped edges contains loads more than the number of kmers in the read, the mappings probably aren't seful
        // actually first just see if it mapping to loads of edges is an issue
        std::vector<edgeKmerPosition> mapped_edges = kMatch.lookupRead(read);
        int distinct_edge_ids = 0;
        if (mapped_edges.size() != 0){
            int current_edge_id = mapped_edges[0].edge_id;
            mapping_counts[current_edge_id] += 1;
            for (auto mapping: mapped_edges){
                if (mapping.edge_id != current_edge_id){
                    current_edge_id = mapping.edge_id;
                    distinct_edge_ids += 1;
                }

            }
            if (distinct_edge_ids == 1){
                mapped_to_single_edge += 1;
            }
            if (distinct_edge_ids > 1){
                mapped_to_multiple_edge += 1;
            }
        }
        this->marker_edge_maps.push_back(mapped_edges);
    }
    std::cout << Date() << "finished read mapping"<< std::endl;
    PrintMemUsage(); // 1231452, so increase in memory usage over the course of this function is 105876
    std::cout << "Unmapped reads r1: " << unmapped_marker << std::endl;
    std::cout << "Mapped to single edge_r1: " << mapped_to_single_edge << std::endl;
    std::cout << "Mapped to multiple edges_r2: " << mapped_to_multiple_edge << std::endl;
    std::cout << "TOtal reads: " << marker_sequences.size() << std:: endl;

}


ReadPath MarkerMapper::getFullyMappedProbes(std::vector<edgeKmerPosition> read_mapping, int read_length, int i, int k=31){
    int kmer_errors_allowed = 0;
    int current_edge_id = read_mapping[0].edge_id;
    int consecutive_offsets = 0;
    int last_offset = read_mapping[0].edge_offset;
    ReadPath paths;
    bool main_graph;
    for (auto it = read_mapping.begin(); it != read_mapping.end(); ++it) {
        main_graph = it->edge_id < inv[it->edge_id] ? true: false;
        // if we're still on the same edge
        if (it->edge_id == current_edge_id) {
            // if the edge offset is consecutive
            if (it->edge_offset == (last_offset + 1)) {
                consecutive_offsets += 1;

            }
        }
            // if we're onto a new edge, then check if we completed the previous edge -
        else {
            if (((std::find(hbv.From(mToRight[current_edge_id]).begin(),hbv.From(mToRight[current_edge_id]).end(), it->edge_id) != hbv.From(mToRight[current_edge_id]).end()) && main_graph) || ((std::find(hbv.From(mToLeft[current_edge_id]).begin(),hbv.From(mToLeft[current_edge_id]).end(), it->edge_id) != hbv.From(mToLeft[current_edge_id]).end())&& !main_graph)){
                current_edge_id = it->edge_id;
            }
            // if the number of consecutive offsets is equal to the number of kmers on the read, -1 because we count transitions its fully mapped, L - k
            else if ((consecutive_offsets == hbv.EdgeObject(current_edge_id).size() - k) || consecutive_offsets == read_length - k ) { // lmp reads are shorter than edges, but check both whether entire read OR entire edge is mapped
                paths.push_back(current_edge_id);
            } else {

                current_edge_id = it->edge_id;
                consecutive_offsets = 0;
            }
        }
        last_offset = it->edge_offset;
        // if the above conditional would miss this because we're at the end of the loop, exit here
        if (std::next(it) == read_mapping.end() && ((consecutive_offsets == hbv.EdgeObject(current_edge_id).size() - k) || consecutive_offsets == read_length - k )){
            paths.push_back(current_edge_id);
        }
    }
    return paths;
}