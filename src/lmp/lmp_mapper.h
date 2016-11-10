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
} LMPPair;

class LMPMapper{
    public:
        LMPMapper(vecbvec lmp_reads, HyperBasevector& hbv, KMatch kmatch);
        vecbvec lmp_reads;
        std::vector<LMPPair > read_paths;
        void LMPMapper::LMPReads2MappedPairedEdgePaths();

    private:
        KMatch kMatch;
        HyperBasevector hbv;
        std::vector<std::vector<edgeKmerPosition> > read_edge_maps;
        std::vector<std::vector<edgeKmerPosition> > initalise_read_edge_map(){std::vector<std::vector<edgeKmerPosition> > res; return res;};
        std::vector<ReadPath> initalise_read_path(){std::vector<ReadPath> res; return res;};
        ReadPath LMPMapper::getFullyMappedEdges(std::vector<edgeKmerPosition> read_mapping, int read_length, int count,  int k=31);
        void LMPMapper::mapReads();
        std::vector<LMPPair > LMPMapper::readEdgeMap2LMPPairs();
         ReadPath LMPMapper::sortMappingsFindFullyMappedEdges(std::vector<edgeKmerPosition>  read_mapping, int read_length, int count);
        void LMPMapper::findFullyMappedEdges();
        bool sanityCheck(LMPPair lmp_pair, int i,  int insert_size=8000);
        void LMPMapper::removeUselessLMPMappings(std::vector<LMPPair > &read_paths, std::vector<LMPPair > &read_paths_for_scaffolding);
        //bool LMPMapper::compareEdgeKmerPositions(const edgeKmerPosition &ekp1, const edgeKmerPosition &ekp2);
};
#endif //W2RAP_CONTIGGER_LMP_MAPPER_H
