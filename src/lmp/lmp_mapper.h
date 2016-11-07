//
// Created by Katie Barr (EI) on 04/11/2016.
//

#ifndef W2RAP_CONTIGGER_LMP_MAPPER_H
#define W2RAP_CONTIGGER_LMP_MAPPER_H

#include <Basevector.h>
#include <kmers/kmatch/KMatch.h>

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
        std::vector<int> LMPMapper::getFullyMappedEdges(std::vector<edgeKmerPosition> read_mapping, int k=31);
        //bool LMPMapper::compareEdgeKmerPositions(const edgeKmerPosition &ekp1, const edgeKmerPosition &ekp2);
};
#endif //W2RAP_CONTIGGER_LMP_MAPPER_H
