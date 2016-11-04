//
// Created by Katie Barr (EI) on 04/11/2016.
//


#include <paths/HyperBasevector.h>
#include "lmp_mapper.h"


LMPMapper::LMPMapper(vecbvec lmp_reads, HyperBasevector& hbv, KMatch kmatch): lmp_reads(lmp_reads), hbv(hbv), kMatch(kmatch) {}

std::vector<std::vector<edgeKmerPosition> > LMPMapper::mapReads(){
    kMatch.Hbv2Map(&hbv);
    std::cout << kMatch.edgeMap.size() << std::endl;
    std::vector<std::vector<edgeKmerPosition> > read_mappings;
    //std::vector<edgeKmerPosition> res = kmatch.lookupRead(lmp_data[0].ToString());
    for (int i=0; i < lmp_reads.size(); i++){
        read_mappings.push_back(kMatch.lookupRead(lmp_reads[i].ToString()));
    }
    return read_mappings;
}