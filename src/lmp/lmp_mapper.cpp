//
// Created by Katie Barr (EI) on 04/11/2016.
//


#include <paths/HyperBasevector.h>
#include <kmers/kmatch/KMatch.h>
#include "lmp_mapper.h"
#include <algorithm>
#include <paths/PathFinder_kb.h>

LMPMapper::LMPMapper(vecbvec* lmp_reads, HyperBasevector& hbv, vec<int>& inv, KMatch kmatch): lmp_reads(lmp_reads), hbv(hbv), inv(inv), kMatch(kmatch), read_edge_maps(initalise_read_edge_map()) {}

bool compareEdgeKmerPositions(const edgeKmerPosition &ekp1, const edgeKmerPosition &ekp2){
    // sort by edge id first, then offset. so a higher edge id with a lower offset would be the greater one
    if (ekp1.edge_id < ekp2.edge_id){
        return true;
    } else if (ekp1.edge_id == ekp2.edge_id && ekp1.edge_offset < ekp2.edge_offset){
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

// stolen from gonza so blame him!!
std::vector<edgeKmerPosition> LMPMapper::readOffsetFilter(std::vector<edgeKmerPosition> data){
    // Count the number of edges_id that map to each kmer in the read, filter out the edges that map to a position that has more than one edge mapped to it
            // input is a vector with all the links for a specific rea
    std::vector<edgeKmerPosition> links;
    std::map<int, unsigned int> pos_count;
    for (auto l: data){
        // Aumentar el contador, valores son inicializados a 0!?
                pos_count[l.read_offset] += 1;
    }
    // Solo seleccionados lo que estan en offsets de la lectura con un solo edge mapeado
    for (auto l: data){
        if (pos_count[l.read_offset] == 1){
                links.push_back(l);
            }
        }
    return links;
}

void LMPMapper::mapReads(){
    kMatch.Hbv2Map(&hbv);
    std::cout << Date() << "starting lmp read mapping"<< std::endl;
    int mapped_to_single_edge_r1 = 0;
    int mappted_to_multiple_edge_r1 = 0;
    int mapped_to_single_edge_r2 = 0;
    int mappted_to_multiple_edge_r2 = 0;
    int unmapped_reads_r1 = 0;
    int unmapped_reads_r2 = 0;
    int mapped_to_single_edge_rc = 0;
    int mappted_to_multiple_edge_rc = 0;
    std::map<uint64_t, int> mapping_counts;
    #pragma omp parallel for
    for (int i=0; i < lmp_reads->size(); i++){
        std::string read = (*lmp_reads)[i].ToString();
        std::vector<edgeKmerPosition> mapped_edges = kMatch.lookupRead(read);
        int distinct_edge_ids = 0;
        if (mapped_edges.size() != 0){
            int distinct_edge_ids = +1;
            int current_edge_id = mapped_edges[0].edge_id;
            mapping_counts[current_edge_id] += 1;
            for (auto mapping: mapped_edges){
                if (mapping.edge_id != current_edge_id){
                    current_edge_id = mapping.edge_id;
                    distinct_edge_ids += 1;
                    mapping_counts[current_edge_id] += 1;

                }
        }
            if (distinct_edge_ids == 1 && i%2==1){
                mapped_to_single_edge_r2 += 1;
            }
            if (distinct_edge_ids > 1 && i%2==1){
                mappted_to_multiple_edge_r2 += 1;
            }
            if (distinct_edge_ids == 1 && i%2==0){
                mapped_to_single_edge_r1 += 1;
            }
            if (distinct_edge_ids > 1 && i%2==0){
                mappted_to_multiple_edge_r1 += 1;
            }
        } else {
            if (mapped_edges.size() == 0 && i%2==0){
            unmapped_reads_r1 += 1;
        } else if (mapped_edges.size() == 0 && i%2==1){
                    unmapped_reads_r2 += 1;

            }
        }
#pragma omp critical (lmp_read_map)
        read_edge_maps.push_back(mapped_edges);
    }
    std::cout << Date() << "finished read mapping"<< std::endl;

    std::cout << "Unmapped reads r1: " << unmapped_reads_r1 << std::endl;
    std::cout << "Unmapped reads r2: " << unmapped_reads_r2 << std::endl;
    std::cout << "Mapped to single edge_r1: " << mapped_to_single_edge_r1 << std::endl;
    std::cout << "Mapped to multiple edges_r1: " << mappted_to_multiple_edge_r1 << std::endl;
    std::cout << "Mapped to single edge_r2: " << mapped_to_single_edge_r2 << std::endl;
    std::cout << "Mapped to multiple edges_r2: " << mappted_to_multiple_edge_r2 << std::endl;
    std::cout << "TOtal reads: " << lmp_reads->size() << std:: endl;
    std::ofstream edge_map_counter;
    edge_map_counter.open("/Users/barrk/Documents/arabidopsis_data/edge_mapping_counts_raw.txt");
    for (auto edge_map_count: mapping_counts){
        edge_map_counter << edge_map_count.first << ', ' << edge_map_count.second << std::endl;
    }
    edge_map_counter.close();

}

void LMPMapper::LMPReads2MappedPairedEdgePaths(std::vector<LMPPair > & lmp_pairs_for_scaffolding, std::vector<LMPPair > & lmp_pairs_for_insert_size_estimation, std::map<uint64_t, std::vector<int> > & edge_id_to_pair_id_map){
    mapReads();
    readEdgeMap2LMPPairs(lmp_pairs_for_scaffolding,  lmp_pairs_for_insert_size_estimation, edge_id_to_pair_id_map);
    std::cout << "LMPReads2MappedPairedEdgePaths edge id map size: " << edge_id_to_pair_id_map.size() << std::endl;

}


std::string LMPMapper::path_str(ReadPath path) {
    std::string s="[";
    for (auto p:path){
        // output edge id on hbv and involution
        s+=std::to_string(p)+":"+std::to_string(inv[p])+" ";//+" ("+std::to_string(mHBV.EdgeObject(p).size())+"bp "+std::to_string(paths_per_kbp(p))+"ppk)  ";
    }
    s+="]";
    return s;
}

void LMPMapper::readEdgeMap2LMPPairs(std::vector<LMPPair >  & lmp_pairs_for_scaffolding, std::vector<LMPPair > & lmp_pairs_for_insert_size_estimation, std::map<uint64_t, std::vector<int> > & edge_id_to_pair_id_map){
    int counter = 0;
    int counter_p1 = 0;
    int counter_p2 = 0;
    int counter_single_reads_mapping_different_edges = 0;
    int counter_both_edges_mapping_to_subgraph_potentially_solveably = 0;
    std::map<std::pair<uint64_t, uint64_t>, int > counts_for_each_solveabe_pair;
    std::map<std::pair<uint64_t, uint64_t>, std::vector<int> > pair_ids_for_each_solveabe_pair;
    std::map<uint64_t, int> mapping_counts;
    int map_to_same_edge=0;
    // to determine if region is solvable, need to
    std::map<uint64_t, std::vector<LMPPair>> edge_id_lmp_dict;
    int pair_id = 0;
    for (int i=0; i < read_edge_maps.size() - 1; ++i) {
        LMPPair lmp_pair;
        lmp_pair.read_index = i / 2;
        std::vector<edgeKmerPosition> read_mapping_p1 = read_edge_maps[i];
        int read_len = static_cast<int>((*lmp_reads)[i].ToString().size());
        // ensure that full edges will be together, with offsets in increasing order, in read_mapping vector
        lmp_pair.p1 = sortMappingsFindFullyMappedEdges(read_mapping_p1, read_len, i);
        i = i + 1;
        std::vector<edgeKmerPosition> read_mapping_p2 = read_edge_maps[i];
        lmp_pair.p2 = sortMappingsFindFullyMappedEdges(read_mapping_p2, read_len, i);
        if (lmp_pair.p1.size() != 0) {
            if (lmp_pair.p2.size() != 0) {

                counter += 1;
                //for (auto edge1: lmp_pair.p1) {
                if (lmp_pair.p1.size() == 1 && lmp_pair.p2.size() == 1) {
                    // pretty sure we do not want to allow mappings where one read is not reverse complemented!!
                    if (lmp_pair.p1[0] == inv[lmp_pair.p2[0]] || lmp_pair.p1[0] == lmp_pair.p2[0]) {
                        // todo: lmp_pair is wrong data type for this, as we need edge offset
                        lmp_pairs_for_insert_size_estimation.push_back(lmp_pair);
                    } else { // if they both map to a single edge, but not the same edge, it can be used to join 2 edges, maybe...

                        lmp_pair.pair_id = pair_id;
                        lmp_pairs_for_scaffolding.push_back(
                                lmp_pair);
                        pair_id += 1;
                        counter_single_reads_mapping_different_edges += 1;
                        for (auto edge1: lmp_pair.p1) {
                            mapping_counts[edge1] += 1;
                        }
                        for (auto edge2: lmp_pair.p2) {
                            mapping_counts[edge2] += 1;
                        }
                        auto edge1 = lmp_pair.p1[0];
                        edge_id_to_pair_id_map[edge1].push_back(lmp_pair.pair_id);
                        auto edge2 = lmp_pair.p2[0];
                        edge_id_to_pair_id_map[edge2].push_back(lmp_pair.pair_id);

                    }
                    counter_p2 += 1;
                }
                counter_p1 += 1;

            }
        }
    }

    std::cout << "LMP reads 1 and 2 map to same egde: " << map_to_same_edge << std::endl;
    std::cout << "total both pairs mapped: " << counter << " p1 mapped: " << counter_p1 << "p2 mapped " << counter_p2 << std::endl;
    std::cout << "readEdgeMap2LMPPairs lmp pairs size: " << lmp_pairs_for_scaffolding.size() << std::endl;
        std::cout << "pairs map to edge and reverse complement and signle edge: "
                  << lmp_pairs_for_insert_size_estimation.size() << std::endl;

}




ReadPath LMPMapper::sortMappingsFindFullyMappedEdges(std::vector<edgeKmerPosition>  read_mapping, int read_length, int i){
    if (read_mapping.size() > 0) {
        std::sort(read_mapping.begin(), read_mapping.end(), compareEdgeKmerPositions);

        return getFullyMappedEdges(read_mapping, read_length, i);
    }
    // prevent bad access when no reads mapped
    ReadPath empty_path;
    return empty_path;
}


// edges are far longer than reads,
// assume perfect mapping to start with
ReadPath LMPMapper::getFullyMappedEdges(std::vector<edgeKmerPosition> read_mapping, int read_length, int i, int k=31){
    bool is_r2 = i % 2 != 0;
    int kmer_errors_allowed = 0;
    if (is_r2){
        // the read 2 error rate at the ends i much higher, so don't consider the last  kmers
        kmer_errors_allowed = (read_length/2) - k;
    }
    int current_edge_id = read_mapping[0].edge_id;
    int consecutive_offsets = 0;
    int last_offset = read_mapping[0].edge_offset;
    ReadPath paths;
    for (auto it = read_mapping.begin(); it != read_mapping.end(); ++it) {
        // if we're still on the same edge
        if (it->edge_id == current_edge_id) {
            // if the edge offset is consecutive
            if (it->edge_offset == (last_offset + 1)) {
                consecutive_offsets += 1;

            }
        }
            // if we're onto a new edge, then check if we completed the previous edge -
        else {
            // if the number of consecutive offsets is equal to the number of kmers on the read, -1 because we count transitions its fully mapped, L - k
            if ((consecutive_offsets == hbv.EdgeObject(current_edge_id).size() - k) || consecutive_offsets == read_length - k || (is_r2 && (consecutive_offsets > read_length - kmer_errors_allowed))) { // lmp reads are shorter than edges, but check both whether entire read OR entire edge is mapped
                paths.push_back(current_edge_id);
            }
            current_edge_id = it->edge_id;
            consecutive_offsets = 0;
        }
        last_offset = it->edge_offset;
        // if the above conditional would miss this because we're at the end of the loop, exit here
        if (std::next(it) == read_mapping.end() && ((consecutive_offsets == hbv.EdgeObject(current_edge_id).size() - k) || consecutive_offsets == read_length - k || (is_r2 && (consecutive_offsets > read_length - k - kmer_errors_allowed)))){
            paths.push_back(current_edge_id);
        }
    }
    return paths;
}


