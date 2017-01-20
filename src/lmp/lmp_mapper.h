//
// Created by Katie Barr (EI) on 04/11/2016.
//

#ifndef W2RAP_CONTIGGER_LMP_MAPPER_H
#define W2RAP_CONTIGGER_LMP_MAPPER_H

#include <Basevector.h>
#include <kmers/kmatch/KMatch.h>
#include <paths/long/ReadPath.h>

// perhaps better to pair all along rather than find pair at the end
typedef struct {
    ReadPath p1;
    ReadPath p2;
    int read_index;
    int pair_id;
} LMPPair;

class LMPMapper{
    public:
        LMPMapper(vecbvec* lmp_reads, HyperBasevector& hbv, vec<int>& inv, KMatch kmatch);
        vecbvec* lmp_reads;
        std::vector<LMPPair > read_paths;
        void LMPReads2MappedPairedEdgePaths(std::vector<LMPPair >  & lmp_pairs_for_scaffolding, std::vector<LMPPair > & lmp_pairs_for_insert_size_estimation, std::map<uint64_t, std::vector<int> > & edge_id_to_pair_id_map);

    private:
        KMatch kMatch;
        HyperBasevector hbv;
        vec<int> inv;
        std::vector<std::vector<edgeKmerPosition> > read_edge_maps;
        std::vector<std::vector<edgeKmerPosition> > initalise_read_edge_map(){std::vector<std::vector<edgeKmerPosition> > res; return res;};
        std::vector<ReadPath> initalise_read_path(){std::vector<ReadPath> res; return res;};
        ReadPath getFullyMappedEdges(std::vector<edgeKmerPosition> read_mapping, int read_length, int count,  int k=31);
        void mapReads();
        void readEdgeMap2LMPPairs(std::vector<LMPPair >  & lmp_pairs_for_scaffolding, std::vector<LMPPair > & lmp_pairs_for_insert_size_estimation, std::map<uint64_t, std::vector<int> > & edge_id_to_pair_id_map);
         ReadPath sortMappingsFindFullyMappedEdges(std::vector<edgeKmerPosition>  read_mapping, int read_length, int count);
        std::vector<edgeKmerPosition> LMPMapper::readOffsetFilter(std::vector<edgeKmerPosition> data);
        std::string path_str(ReadPath path);
        };
#endif //W2RAP_CONTIGGER_LMP_MAPPER_H
