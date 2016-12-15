//
// Created by Katie Barr (EI) on 15/12/2016.
//

#include "ComplexRegion.h"
using namespace LMPRegions;


ComplexRegion::ComplexRegion(){};

ComplexRegion::ComplexRegion(std::vector<uint64_t  > edges_in, std::vector<uint64_t  > edges_out, vec<int>& involution, int insert_size=5000):
    edges_in(edges_in), edges_out(edges_out), insert_size(insert_size), involution(involution)
{};


void ComplexRegion::AddPath(ReadPath path){
    std::vector<uint64_t> path_canonical = canonicalisePath(path);
    // sanity check path, ensure its not a cycle, and is in this region
    if (path.front() != path.back() &&
            std::find(edges_in_canonical.begin(), edges_in_canonical.end(), path_canonical.front()) != edges_in_canonical.end()
            && std::find(edges_out_canonical.begin(), edges_out_canonical.end(), path_canonical.back()) != edges_out_canonical.end()){
        candidate_paths.push_back(path_canonical);
    } else {
        std::cout << "Path cannot be added to region" << std::endl;
    }
}


std::vector<uint64_t>  ComplexRegion::canonicalisePath(ReadPath path){
    std::vector<uint64_t> result;
    bool involution_path = false;
    // if first edge is on involution, all will be
    if (path[0] > involution[path[0]]){
        involution_path = true;
    }
    for (auto edge: path){
        if (involution_path){
            result.push_back(involution[edge]);
        } else {
            result.push_back(edge);
        }
    }
    if (involution_path || result.back() < result.front()) {
        std::reverse(result.begin(), result.end());
    }
    return result;

}
