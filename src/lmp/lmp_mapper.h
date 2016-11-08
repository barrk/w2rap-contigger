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
    ReadPath p1; // TODO ask gonza more informative name for this
    ReadPath p2; // TODO ask gonza more informative name for this
} LMPPair;

class LMPMapper{
    public:
        LMPMapper(vecbvec lmp_reads, HyperBasevector& hbv, KMatch kmatch);
        vecbvec lmp_reads;
        void LMPMapper::mapReads();
        void LMPMapper::findFullyMappedEdges();

    private:
        KMatch kMatch;
        HyperBasevector hbv;
        std::vector<std::vector<edgeKmerPosition> > read_edge_maps;
        std::vector<std::vector<edgeKmerPosition> > initalise_read_edge_map(){std::vector<std::vector<edgeKmerPosition> > res; return res;};
        std::vector<LMPPair > read_paths;
        std::vector<ReadPath> initalise_read_path(){std::vector<ReadPath> res; return res;};
        ReadPath LMPMapper::getFullyMappedEdges(std::vector<edgeKmerPosition> read_mapping, int k=31);
        void LMPMapper::convertLMPPairsToReadPaths();
         ReadPath LMPMapper::sortMappingsFindullyMappedEdges(std::vector<edgeKmerPosition>);
        // same edge wooll occur in multiple paths- hot to determine best choice?
        //ReadPath LMPMapper::getReadMathPair(int edge_id);
        //ReadPath LMPMapper::getReadMathPair(int read_index);
        //bool LMPMapper::compareEdgeKmerPositions(const edgeKmerPosition &ekp1, const edgeKmerPosition &ekp2);
};
#endif //W2RAP_CONTIGGER_LMP_MAPPER_H
