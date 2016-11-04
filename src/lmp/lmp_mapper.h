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
        std::vector<std::vector<edgeKmerPosition> > LMPMapper::mapReads();


    private:
        KMatch kMatch;
        HyperBasevector hbv;
};
#endif //W2RAP_CONTIGGER_LMP_MAPPER_H
