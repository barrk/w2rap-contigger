//
// Created by Katie Barr (EI) and Rob Vickerstaff (EMR) on 24/02/2017.
//

#ifndef W2RAP_CONTIGGER_MARKER_MAPPING_H
#define W2RAP_CONTIGGER_MARKER_MAPPING_H

#include <Basevector.h>
#include <kmers/kmatch/KMatch.h>
#include <paths/long/ReadPath.h>

class MarkerMapper{
public:
    MarkerMapper(vecbvec& lmp_reads, std::vector<std::string> marker_ids, HyperBasevector& hbv, vec<int>& inv, KMatch kmatch);
    vecbvec& marker_sequences;

private:
    KMatch kMatch;
    HyperBasevector hbv;
    vec<int> inv;
    vec<int> mToLeft;
    vec<int> mToRight;
    std::vector<std::string> marker_ids;
    std::vector<std::vector<edgeKmerPosition> > marker_edge_maps;
    std::vector<std::vector<edgeKmerPosition> > initalise_marker_edge_map(){std::vector<std::vector<edgeKmerPosition> > res; return res;};
    void mapMarkers();
    ReadPath getFullyMappedProbes(std::vector<edgeKmerPosition> marker_mapping, int read_length, int i, int k=31);

};
#endif //W2RAP_CONTIGGER_MARKER_MAPPING_H
