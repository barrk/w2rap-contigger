#include "KMatch.h"
#include <sys/time.h>
#include <thread>
#include <paths/HyperBasevector.h>

KMatch::KMatch(int kv){
  if (kv > 31){
    std::cout << "Kmer value is too big for this, using 31 instead... " << std::endl;
    this->K = 31;
  }
  this->K = kv;
}

std::vector<pKmer> KMatch::ProduceKmers(const std::string seq) const {
  // get a sequence a produce the set of kmers ()
  std::vector<pKmer> kmer_vector;

  int offset=0;
  const uint64_t KMER_MASK=( ((uint64_t)1)<<(this->K*2) )-1;

  if (seq.size()>this->K) {
    const char *s = seq.c_str();
    int64_t last_unknown = -1;
    uint64_t fkmer = 0;
    for (auto p=0; p < seq.size(); ++p) {
      //fkmer: grows from the right (LSB)
      switch (s[p]) {
        case 'A':
        case 'a':
          fkmer = ((fkmer << 2) + 0) & KMER_MASK;
          break;
        case 'C':
        case 'c':
          fkmer = ((fkmer << 2) + 1) & KMER_MASK;
          break;
        case 'G':
        case 'g':
          fkmer = ((fkmer << 2) + 2) & KMER_MASK;
          break;
        case 'T':
        case 't':
          fkmer = ((fkmer << 2) + 3) & KMER_MASK;
          break;
        default:
          fkmer = ((fkmer << 2) + 0) & KMER_MASK;
          last_unknown = p;
          break;
      }
      //TODO: last unknown passed by?
      if (last_unknown + this->K <= p) {
        pKmer pair_temp;
        pair_temp.kmer = fkmer;
        pair_temp.offset = offset;
        kmer_vector.push_back(pair_temp);
        offset++;
      }
    }
  }
  return kmer_vector;
}

void KMatch::Hbv2Map(HyperBasevector* hbv){
  //  std::vector<kmer_position_t> karray;

  std::map<uint64_t, std::vector<edgeKmerPosition>> edgeDict;
  uint32_t seq_index=0;

  auto edges = hbv->Edges();
  std::vector<edgeKmerPosition> temp_vector;
  edgeKmerPosition tmatch;
  for (auto seqN=0; seqN<edges.size(); ++seqN) {
    auto seq = edges[seqN].ToString();
    auto kv = this->ProduceKmers(seq);

    for (auto a=0; a<kv.size(); ++a){
      temp_vector.clear();
      if (this->edgeMap.find(kv[a].kmer) == this->edgeMap.end()){

        tmatch.edge_id = seq_index;
        tmatch.edge_offset = kv[a].offset;
        temp_vector.push_back(tmatch);
        this->edgeMap[kv[a].kmer] = temp_vector;

      } else {
        auto temp_vector = this->edgeMap[kv[a].kmer];
        tmatch.edge_id = seq_index;
        tmatch.edge_offset = kv[a].offset;
        temp_vector.push_back(tmatch);
        this->edgeMap[kv[a].kmer] = temp_vector;
      }

    }
    seq_index++;
  }
  std::cout << Date() << "finished building edge dict"<< std::endl;

}

std::vector<edgeKmerPosition> KMatch::lookupRead(const std::string read) const {
  // produce kmers
  auto read_kmers = this->ProduceKmers(read);
    //std::cout << "mapping read: " << read << std::endl;
  // look kmers in the dictionary
  std::vector<edgeKmerPosition> mapped_edges;
  int cont = 0; // Cont to hold the kmer offset in the read
  for (auto const &kmer: read_kmers){
    std::map<uint64_t, std::vector<edgeKmerPosition>>::const_iterator matched_edges = this->edgeMap.find(kmer.kmer);
    if (matched_edges != this->edgeMap.end()){
      for (auto edge_position: matched_edges->second){
        edgeKmerPosition x;
        x.edge_id = edge_position.edge_id;
        x.edge_offset = edge_position.edge_offset;
        x.read_offset = cont;
        x.kmer = edge_position.kmer;
        mapped_edges.push_back(x);
      }
    }
    cont++;
  }
  return mapped_edges;
}
