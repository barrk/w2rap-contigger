//
// Created by Katie Barr (EI) and Rob Vickerstaff (EMR) on 24/02/2017.
//

#include "marker_mapping.h"

MarkerMapper::MarkerMapper(vecbvec& marker_sequences, HyperBasevector& hbv, vec<int>& inv, KMatch kmatch): marker_sequences(marker_sequences), hbv(hbv), kMatch(kmatch) {}



void MarkerMapper::mapMarkers(){
    kMatch.Hbv2Map(&hbv);
    PrintMemUsage(); // 1125576
    std::cout << Date() << "starting lmp read mapping"<< std::endl;
    int mapped_to_single_edge = 0;
    int mapped_to_multiple_edge = 0;
    int unmapped_marker = 0;
    std::map<uint64_t, std::atomic<int>> mapping_counts;
    for (int i=0; i < marker_sequences.size(); i++){
        std::string read = (marker_sequences)[i].ToString();
        // make guess that if mapped edges contains loads more than the number of kmers in the read, the mappings probably aren't seful
        // actually first just see if it mapping to loads of edges is an issue
        std::vector<edgeKmerPosition> mapped_edges = kMatch.lookupRead(read);
        if (mapped_edges.size() > 1000){
            std::cout << "mapped edges size: " << mapped_edges.size() << " read index: " << i << std::endl;
        }
        int distinct_edge_ids = 0;
        if (mapped_edges.size() != 0){
            int current_edge_id = mapped_edges[0].edge_id;
            mapping_counts[current_edge_id] += 1;
            for (auto mapping: mapped_edges){
                if (mapping.edge_id != current_edge_id){
                    current_edge_id = mapping.edge_id;
                    distinct_edge_ids += 1;
                }
                if (distinct_edge_ids > 5){// total guess....
                    mapped_edges.clear();
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