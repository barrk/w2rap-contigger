//
// Created by Katie Barr (EI) on 15/12/2016.
//

#include "ComplexRegion.h"


ComplexRegion::ComplexRegion(){};

ComplexRegion::ComplexRegion(std::vector<uint64_t  > edges_in, std::vector<uint64_t  > edges_out, vec<int>& involution, int insert_size=5000):
    edges_in(edges_in), edges_out(edges_out), insert_size(insert_size), involution(involution)
{
    std::sort(edges_in.begin(), edges_in.end());
    std::sort(edges_out.begin(), edges_out.end());
    edges_in_detailed.resize(edges_in.size());
    edges_out_detailed.resize(edges_out.size());
    for (int i=0; i < edges_in.size(); i++){
        edges_in_detailed[i].edge_id = edges_in[i];
    }
    for (int i=0; i < edges_out.size(); i++) {
        edges_out_detailed[i].edge_id = edges_out[i];
    }
    std::cout << "Adding region with edges in: " << path_str(edges_in) << std::endl;
    std::cout << "Adding region with edges out: " << path_str(edges_out) << std::endl;
    solved = false;
};



std::string ComplexRegion::path_str(std::vector<uint64_t> path) {
    std::string s="[";
    for (auto p:path){
        // output edge id on hbv and involution
        s+=std::to_string(p)+":"+std::to_string(involution[p])+" ";
    }
    s+="]";
    return s;
}

void ComplexRegion::AddPairId(int edge_index, int pair_id, bool in){
    BoundingEdge e;
    if (in){
        e = edges_in_detailed[edge_index];
    } else{
        e = edges_out_detailed[edge_index];

    }
    e.pair_ids.push_back(pair_id);
    if (in){
        edges_in_detailed[edge_index] = e;
    } else{
        edges_out_detailed[edge_index] = e;

    }
}

void ComplexRegion::AddPathTo(int edge_index, bool in, std::vector<uint64_t> path_from_center){
    BoundingEdge e;
    if (in){
        e = edges_in_detailed[edge_index];
    } else{
        e = edges_out_detailed[edge_index];

    }
    e.path_from_center = path_from_center;
    if (in){
        edges_in_detailed[edge_index] = e;
    } else{
        edges_out_detailed[edge_index] = e;

    }
}

void ComplexRegion::FindSolveablePairs(){
    for (auto edge_in: edges_in_detailed) {
        auto in_ids = edge_in.pair_ids;
        for (auto in_id: in_ids){
            for (auto edge_out: edges_out_detailed) {
                auto out_ids = edge_out.pair_ids;
                if ((std::find(out_ids.begin(), out_ids.end(), in_id) !=
                        out_ids.end())) {
                    auto key = std::make_pair(edge_in.edge_id, edge_out.edge_id);
                    combination_counts[key] += 1;

                }
            }
        }
    }
}





void ComplexRegion::isSolved(int min_count){
    // this assumes there is exactly one path per in/out combination- is this always the case? if there are fewer, region is not solved, if there are more over min count, region is not unambiguously solved, though in practise if one in/out pair had much more support than the other, we'd want that
    //check for duplication - i orginally used a set instead of a vector, but changed it for a reson i now can't remember
    std::set<uint64_t> in_edges_solved_set;
    std::set<uint64_t> out_edges_solved_set;
    for (auto count:combination_counts){
        if (count.second > min_count){
            auto in = count.first.first;
            auto out = count.first.second;
            in_edges_solved.push_back(in);
            out_edges_solved.push_back(out); //TODO: think about removing duplication of data
            in_edges_solved_set.insert(in);
            out_edges_solved_set.insert(out);
        }
    }

    std::cout << "in edges solved size: " << in_edges_solved.size() << "edges in size: " <<  edges_in.size() <<  "in edges solved set size: " << in_edges_solved_set.size() << std::endl;
    std::cout << "out edges solved size: " << out_edges_solved.size() << "edges out size: " <<  edges_out.size() <<  "out edges solved set size: " << out_edges_solved_set.size() << std::endl;

    // think about if this is actually the case....
    if((in_edges_solved.size() == edges_in.size()) && (in_edges_solved.size()  == in_edges_solved_set.size()) && (out_edges_solved.size() == edges_out.size()) && (out_edges_solved.size()  == out_edges_solved_set.size())){
    //if(in_edges_solved.size() == edges_in.size()  && out_edges_solved.size() == edges_out.size() ){

        for (int i = 0; i < in_edges_solved.size(); i++){
            combinations_to_use.push_back(std::make_pair(in_edges_solved[i], out_edges_solved[i]));
        }
        solved = true;
        std::cout << "Combinations to use size of region for solved region: " << combinations_to_use.size() << std::endl;
        std::cout << "Region solved" << std::endl;
    }
    std::cout << "Combinations to use size of region: " << combinations_to_use.size() << std::endl;
}

std::vector<std::vector<uint64_t> >  ComplexRegion::BuildPaths(){
    /*
     * determine which in/out pairs to use and build th associated paths
    */
    std::vector<std::vector<uint64_t> > paths;
    std::cout << "building for region, edges in: " << path_str(edges_in) << std::endl;
    std::cout << "building for region, edges out: " << path_str(edges_out) << std::endl;
    // combination counts contains pairs of in/out edge ids,
    for (auto in_out_pair:combinations_to_use){
        auto in_edge_id = in_out_pair.first;
        auto out_edge_id = in_out_pair.second;
        std::cout << "Building path from: " << in_edge_id << " to: " << out_edge_id << " paths size: " << paths.size() << std::endl;
        for (auto in_edge: edges_in_detailed){
            if (in_edge.edge_id == in_edge_id) {
                for (auto out_edge: edges_out_detailed) {
                    if (out_edge.edge_id == out_edge_id){
                        auto path = BuildPath(in_edge, out_edge);
                        paths.push_back(path);
                    }
                }
            }
        }
        std::cout << " paths size: " << paths.size() << std::endl;
    }
    return paths;
}

// this is baed on paths to original edges in/out
std::vector<uint64_t>  ComplexRegion::BuildPath(BoundingEdge edge_in, BoundingEdge edge_out){
    std::vector<uint64_t> result;
    result.push_back(edge_in.edge_id); // should be same as edge id- if this is on inv, all others should be?
    for (auto e:edge_in.path_from_center){
        result.push_back(e);
        std::cout << e << std::endl;
    }
    for (auto e:edge_out.path_from_center){
        result.push_back(e);
        std::cout << e << std::endl;
    }
    result.push_back(edge_out.edge_id);
    return result;

}


ComplexRegionCollection::ComplexRegionCollection(vec<int>& involution): involution(involution){}



bool ComplexRegionCollection::AddRegion(std::vector<uint64_t> edges_in, std::vector<uint64_t> edges_out,
               vec<int> &involution, int insert_size = 5000){
    std::set<uint64_t > check_edges_distinct;
    for (auto edge: edges_in){
        check_edges_distinct.insert(edge);
        check_edges_distinct.insert(involution[edge]);
    }
    for (auto edge: edges_out){
        check_edges_distinct.insert(edge);
        check_edges_distinct.insert(involution[edge]);
    }
    if (check_edges_distinct.size() == (2*edges_in.size() + 2*edges_out.size())) {
        ComplexRegion complex_region(edges_in, edges_out, involution, insert_size);
        complex_regions.push_back(complex_region);
        auto key = std::make_pair(complex_region.edges_in, complex_region.edges_out);
        edges_to_region_index[key] = complex_regions.size();
        return true;
    } else {
        std::cout << "In/Out edges do not define a valid region" << std::endl;
        return false;
    }
}

bool ComplexRegionCollection::ContainsRegionWithEdges(std::vector<uint64_t> edges_in, std::vector<uint64_t> edges_out){
    return edges_to_region_index.count(std::make_pair(edges_in, edges_out)) != 0;
}

std::pair<ComplexRegion, int> ComplexRegionCollection::GetRegionWithEdges(std::vector<uint64_t> edges_in, std::vector<uint64_t> edges_out){
    auto index = edges_to_region_index[std::make_pair(edges_in, edges_out)];
    return std::make_pair(complex_regions[index], index);

}


void ComplexRegionCollection::SelectRegionsForPathSeparation(){
    for (auto region:complex_regions){
        if (region.solved){
            solved_regions.push_back(region);
            std::cout << "Combinations to use size of solved region: " << region.combinations_to_use.size() << std::endl;
        }
    }
    for (int i = 0; i < solved_regions.size() ; i++){
        auto region = solved_regions[i];
        bool region_clashed = false;
        auto edges_in = region.edges_in;
        for (int j = i + 1; j < solved_regions.size(); j++){
            auto region_2 = solved_regions[j];
            for (auto edge: edges_in) {
                // check if it doesn't share an edge with another region, or that the involution, so in the opposite direction,
                if (std::find(region_2.edges_in.begin(), region_2.edges_in.end(), edge) != region_2.edges_in.end()
                        || std::find(region_2.edges_out.begin(), region_2.edges_out.end(), involution[edge]) != region_2.edges_out.end()){
                        solved_regions_final.push_back(FindBestSolvedRegion(region, region_2));
                        region_clashed = true;
                }
            }
        }
        if (!region_clashed){
            solved_regions_final.push_back(region);
        }
    }
    if (solved_regions_final.size() > 1) {
        for (int i = 0; i < solved_regions_final.size(); i++) {
            bool region_removed = false;
            auto region = solved_regions_final[i];
            auto edges_out = region.edges_out;
            for (int j = i + 1; j < solved_regions_final.size() -1; j++) {
                auto region_2 = solved_regions_final[j];
                for (auto edge: edges_out) {
                    if (std::find(region_2.edges_out.begin(), region_2.edges_out.end(), edge) != region_2.edges_out.end()
                        || std::find(region_2.edges_out.begin(), region_2.edges_out.end(), involution[edge]) !=
                           region_2.edges_out.end()) {
                        if (region.edges_in.size() > region_2.edges_in.size()){ //TODO define proper = operator so we cna ue the select best region function here tooo
                            solved_regions_final.erase(solved_regions_final.begin() + j);
                        } else {
                            solved_regions_final.erase(solved_regions_final.begin() + i);
                            region_removed = true; // start from next region if we've removed the region we are currently looping over
                            continue;
                        }
                    }
                }
                if (region_removed){
                    continue;
                }
            }

            if (region_removed){
                continue;
            }
        }
    }
    std::cout << "Number of potential solved regions after removing clashing ones: " << solved_regions_final.size() << std::endl;
}


std::vector<std::vector<uint64_t> > ComplexRegionCollection::GetPathsToSeparate(){
    // for each solved region, get each path in the right direction (so from in to out), on the main graph, not the involution
    std::vector<std::vector<uint64_t> > paths_to_separate;
    std::cout << "number of solved regions: " << solved_regions_final.size()  << std::endl;
    for (auto region:solved_regions_final){
        auto paths = region.BuildPaths();
        for (auto path:paths){
            paths_to_separate.push_back(path);
            std::cout << "adding path: " << path_str(path) <<std::endl;
        }
    }
    std::cout << "patsh to separate, about to return: " << paths_to_separate.size() << " paths." << std::endl;
    return paths_to_separate;
}



std::string ComplexRegionCollection::path_str(std::vector<uint64_t> path) {
    std::string s="[";
    for (auto p:path){
        // output edge id on hbv and involution
        s+=std::to_string(p)+":"+std::to_string(involution[p])+" ";//+" ("+std::to_string(mHBV.EdgeObject(p).size())+"bp "+std::to_string(paths_per_kbp(p))+"ppk)  ";
    }
    s+="]";
    return s;
}

// select one which solves most paths
ComplexRegion ComplexRegionCollection::FindBestSolvedRegion(ComplexRegion region_1, ComplexRegion region_2){
    if (region_1.edges_in.size() > region_2.edges_in.size()){
        return region_1;
    } else {
        return region_2;
    }
}
