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

//TODO repeat all this on involution graph

void LMPMapper::mapReads(){
    kMatch.Hbv2Map(&hbv);
    std::cout << kMatch.edgeMap.size() << std::endl;
    //std::vector<edgeKmerPosition> res = kmatch.lookupRead(lmp_data[0].ToString());
    for (int i=0; i < lmp_reads.size(); i++){
        std::vector<edgeKmerPosition> mapped_edges = kMatch.lookupRead(lmp_reads[i].ToString());
        //std::cout << "Mapped read:" << i << " string " << lmp_reads[i].ToString() << " to " << mapped_edges.size() << "edges" << std::endl;
        read_edge_maps.push_back(mapped_edges);
    }
}

void LMPMapper::LMPReads2MappedPairedEdgePaths(std::vector<LMPPair > & lmp_pairs_for_scaffolding){
    mapReads();
    readEdgeMap2LMPPairs(lmp_pairs_for_scaffolding);
    std::cout << "LMPReads2MappedPairedEdgePaths lmp pairs size: " << lmp_pairs_for_scaffolding.size() << std::endl;

}

void LMPMapper::readEdgeMap2LMPPairs(std::vector<LMPPair > & lmp_pairs_for_scaffolding){
    std::vector<LMPPair > read_paths;
    int counter = 0;
    int counter_p1 = 0;
    int counter_p2 = 0;
    for (int i=0; i < read_edge_maps.size() - 1; ++i){
        LMPPair lmp_pair;
        std::vector<edgeKmerPosition> read_mapping_p1 = read_edge_maps[i];
        int read_len = static_cast<int>(lmp_reads[i].ToString().size());
        // ensure that full edges will be together, with offsets in increasing order, in read_mapping vector
        lmp_pair.p1 = sortMappingsFindFullyMappedEdges(read_mapping_p1, read_len, i);
        i = i + 1;
        std::vector<edgeKmerPosition> read_mapping_p2 = read_edge_maps[i];
        lmp_pair.p2 = sortMappingsFindFullyMappedEdges(read_mapping_p2, read_len, i);
        if (lmp_pair.p1.size() != 0){
            if (lmp_pair.p2.size() != 0){
                counter += 1;
            }
            counter_p1 += 1;
        }
        if (lmp_pair.p2.size() != 0){
            counter_p2 += 1;
        }
        if (!sanityCheck(lmp_pair, i)){
            // decide what to do in this case... sanity check is more useful in working out mapping hasn't gone horriblywrong
            //std::cout << "pair failed sanity check" << std::endl;
            // for now just discard these read paths
        } else {
            read_paths.push_back(
                    lmp_pair);
        }
    }
    std::cout << "total both pairs mapped: " << counter << " p1 mapped: " << counter_p1 << "p2 mapped " << counter_p2 << std::endl;
    removeUselessLMPMappings(read_paths, lmp_pairs_for_scaffolding);
    std::cout << "readEdgeMap2LMPPairs lmp pairs size: " << lmp_pairs_for_scaffolding.size() << std::endl;

}


void LMPMapper::removeUselessLMPMappings(std::vector<LMPPair > &read_paths, std::vector<LMPPair > &read_paths_for_scaffolding){
    // if mappings are both to the same edge, or there are no mappings, remove these from consideration
    for (auto read_path_pair: read_paths){
        if (read_path_pair.p1.size() != 0 && read_path_pair.p2.size() != 0){
            //std::cout << "read path pair has nonzero mappings:" << read_path_pair.p1[0] << read_path_pair.p2[0] << std::endl;
            if (read_path_pair.p1 != read_path_pair.p2){
                read_paths_for_scaffolding.push_back(read_path_pair);
            }

        }

    }
    std::cout << read_paths_for_scaffolding.size() << " of " << read_paths.size() << "lmp pair paths useable for scaffolding" << std::endl;

}



bool LMPMapper::sanityCheck(LMPPair lmp_pair, int i, int insert_size=8000){
    // check that if reads map to long edges, then the pairs map to the same edge, hardcode 'long' to longer than any lmp insert size
    // think conditionals are evaluated left to right, so this will sotp evaluating as soon as one of these is false
   // if both reads map to one edge, and that is a long edge
    if ((lmp_pair.p1.size() == 1 && lmp_pair.p2.size() == 1) && (hbv.EdgeObject(lmp_pair.p1[0]).size() > (insert_size*2))){
        //edgeKmerPosition e1 = read_edge_maps[i-1];
        //edgeKmerPosition e2 = read_edge_maps[i];
        //int edge_sequence_length = static_cast<int>(hbv.EdgeObject(lmp_pair.p1[0]).size());
        // if the read maps to within $insert_size of edge of edge, we don't expect both reads to map to same edge
        //bool offsets_at_edge_of_edge = (e1.offset < (edge_sequence_length - insert_size) || e2.offset < (edge_sequence_length - insert_size) || e1.offset < insert_size || e2.offset < insert_size);
        //if (lmp_pair.p1[0] == lmp_pair.p2[0] || offsets_at_edge_of_edge)){
            return true;
        //}
    } else if (lmp_pair.p1.size() > 1 || lmp_pair.p2.size() > 1 || lmp_pair.p1.size() == 0 || lmp_pair.p2.size() == 0){
        return true;
    }
    return false;
}


ReadPath LMPMapper::sortMappingsFindFullyMappedEdges(std::vector<edgeKmerPosition>  read_mapping, int read_length, int i){
    if (read_mapping.size() > 0) {
        std::sort(read_mapping.begin(), read_mapping.end(), compareEdgeKmerPositions);
        //std::cout << "Read mapping edge id:" << read_mapping[0].edge_id << " offset " << read_mapping[0].offset << " size: " << hbv.EdgeObject(read_mapping[0].edge_id).ToString().length() << " to read of length " << read_length
//                  << std::endl;
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
        // the read 2 error rate at the ends i much higher, so don't consider the last 25 kmers
        kmer_errors_allowed = (read_length/2)-k;
    }
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
                /*if (is_r2 && consecutive_offsets > 60) {
                    std::cout << "edge offset: " << last_offset << " consecutiv offsets: " << consecutive_offsets << " kmer errors allowed " << kmer_errors_allowed << " threshold:  " << read_length - k - kmer_errors_allowed << std::endl;
                    std::cout << "Mapping edge: " << hbv.EdgeObject(current_edge_id).ToString().substr(last_offset, 150) << std::endl;
                    std::cout << "to read length" << read_length << " offset " << last_offset << " string " << lmp_reads[i].ToString() << std::endl;
                }*/
                //std::cout << "incrementing consecutive offsets:" << std::endl;

            }
        }
            // if we're onto a new edge, then check if we completed the previous edge -
        else {
           // std::cout << "current consecutively mapped edges: " << consecutive_offsets << "mapping edge: " << hbv.EdgeObject(it->edge_id).ToString() << std::endl;
            // if the number of consecutive offsets is equal to the number of kmers on the read, -1 because we count transitions its fully mapped, L - k
            if ((consecutive_offsets == hbv.EdgeObject(current_edge_id).size() - k) || consecutive_offsets == read_length - k || (is_r2 && (consecutive_offsets > read_length - k - kmer_errors_allowed))) { // lmp reads are shorter than edges, but check both whether entire read OR entire edge is mapped
                //std::cout << "edge mapped:" << it->edge_id << std::endl;
                paths.push_back(current_edge_id);
            }
            current_edge_id = it->edge_id;
            consecutive_offsets = 0;
        }
        last_offset = it->offset;
        // if the above conditional would miss this because we're at the end of the loop, exit here
        if (std::next(it) == read_mapping.end() && ((consecutive_offsets == hbv.EdgeObject(current_edge_id).size() - k) || consecutive_offsets == read_length - k || (is_r2 && (consecutive_offsets > read_length - k - kmer_errors_allowed)))){
            //std::cout << "edge mapped at end:" << it->edge_id << std::endl;
            paths.push_back(current_edge_id);
        }
    }
    //std::cout << "length of path returned: " << paths.size() << std::endl;
    return paths;
}


