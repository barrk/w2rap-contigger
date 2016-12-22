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
    edges_in_canonical.resize(edges_in.size());
    edges_out_canonical.resize(edges_out.size());
    std::cout << "edges in detailed size constructor: " << edges_in_detailed.size() << std::endl;
    std::cout << "edges out detailed size constructor: " << edges_out_detailed.size() << std::endl;
    std::map<uint64_t, std::vector<int> > pair_ids;
    std::map<uint64_t, uint64_t>  edge_translations;
    canonicaliseEdgesInOut();

};

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
    std::cout << "Adding path: " << std::endl;
    for (auto edge: path_from_center){
        std::cout << " edge " << edge;
    }
    std::cout << std::endl;
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
    std::cout << "finding solveable pairs" << std::endl;
    for (auto edge_in: edges_in_detailed) {
        //std::cout << "Edge in: " << edge_in.edge_id << "pair_ids size: " << edge_in.pair_ids.size() << std::endl;
        auto in_ids = edge_in.pair_ids;
        for (auto in_id: in_ids){
            for (auto edge_out: edges_out_detailed) {
                //std::cout << "Edge out: " << edge_out.edge_id << " in id: " << in_id << std::endl;
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



void  ComplexRegion::canonicaliseEdgesInOut(){
    // this needs to be based on direction of flow through graph
    // to avoid confusion and errors due to reverse complements, and order of read mapping, always deal with paths on the same strand, in the same direction
    //first sort in/out edges- actually maybe should do this at the end
    std::pair<std::vector<uint64_t>, std::vector<BoundingEdge> > in = CanonicaliseEdgeList(edges_in, edges_in_canonical, edges_in_detailed);
    edges_in_canonical = in.first;
    edges_in_detailed = in.second;
    std::pair<std::vector<uint64_t>, std::vector<BoundingEdge> > out = CanonicaliseEdgeList(edges_out, edges_out_canonical, edges_out_detailed);
    edges_out_canonical = out.first;
    edges_out_detailed = out.second;
}

std::pair<std::vector<uint64_t>, std::vector<BoundingEdge> > ComplexRegion::CanonicaliseEdgeList(std::vector<uint64_t> edges, std::vector<uint64_t> edges_canonical, std::vector<BoundingEdge>  detailed_edge_list){
    std::sort(edges.begin(), edges.end());
    for (int i = 0; i < edges.size(); i++){
        auto edge = edges[i];
        BoundingEdge edge_details = detailed_edge_list[i];
        // in practise i think all edges in will be edges in canonical if one of them is, same for out
        if (edge < involution[edge]){// nb edge in i and edge out i may not be together finally, but the purpose of this is to ensure internal consistency
            //std::cout << "edge:" << edge << " is canonical, involution: " << involution[edge] << std::endl;
            edge_details.edge_id = edge;
            edge_details.forward = true;
            edge_details.translated_edge_id = edge;
            edges_canonical[i] = edge;
        } else {
            //std::cout << "edge:" << edge << " canonical version " << edge<< std::endl;
            edge_details.edge_id = involution[edge];
            edge_details.forward = false;
            edge_details.translated_edge_id = involution[edge];
            edges_canonical[i] = involution[edge];

        }
        detailed_edge_list[i] = edge_details;
    }
    return std::make_pair(edges_canonical, detailed_edge_list);
}


void ComplexRegion::isSolved(int min_count){
    // this assumes there is exactly one path per in/out combination- is this always the case? if there are fewer, region is not solved, if there are more over min count, region is not unambiguously solved, though in practise if one in/out pair had much more support than the other, we'd want that
    int number_solved_pairs = 0;
    for (auto count:combination_counts){
        std::cout << "combination: " << count.first.first << ", " << count.first.second << std::endl;
        if (count.second > min_count){
            auto in = count.first.first;
            auto out = count.first.second;
            in_edges_solved.push_back(in);
            out_edges_solved.push_back(out); //TODO: think about removing duplication of data
            combinations_to_use.push_back(std::make_pair(in, out));
            number_solved_pairs += 1;
            std::cout << "solved" << std::endl;
        }
    }
    // think about if this is actually the case....
    if(in_edges_solved.size() == edges_in.size() && out_edges_solved.size() == edges_out.size()){
        solved = true;
        std::cout << "Region solved" << std::endl;
    }
    std::cout << "Number of elements in combination counts: " << combination_counts.size() << std::endl;
}

std::vector<std::vector<uint64_t> >  ComplexRegion::BuildPaths(){
    /*
     * determine which in/out pairs to use and build th associated paths
    */
    std::vector<std::vector<uint64_t> > paths;
    // in edges solved also has nothing in it by here!
    std::cout << "Number of elements in combination counts, build paths: " << combination_counts.size() << std::endl;
    int in_counter = 0;
    int out_counter = 0;
    // combination counts contains pairs of in/out edge ids,
    for (auto in_out_pair:combinations_to_use){
        auto in_edge_id = in_out_pair.first;
        auto out_edge_id = in_out_pair.second;
        // doesn't print this
        std::cout << "Looing for path between " << in_edge_id << " and " << out_edge_id << std::endl;
        for (auto in_edge: edges_in_detailed){
            if (in_edge.edge_id == in_edge_id) {
                std::cout << "found in edge " << std::endl;
                for (auto out_edge: edges_out_detailed) {
                    if (out_edge.edge_id == out_edge_id){
                        auto path = BuildPath(in_edge, out_edge, in_counter, out_counter);
                        paths.push_back(path);
                    }
                    out_counter += 1;
                }
                out_counter = 0;
            }
            in_counter += 1;
        }
        in_counter = 0;

    }
    return paths;
}
// this is baed on paths to original edges in/out
std::vector<uint64_t>  ComplexRegion::BuildPath(BoundingEdge edge_in, BoundingEdge edge_out, int in_counter, int out_counter){
    std::cout << "builng path between: "<< edge_in.edge_id <<" and " <<edge_out.edge_id << std::endl;
    std::vector<uint64_t> result;
    if (edge_in.forward){
        result.push_back(edge_in.translated_edge_id); // should be same as edge id- if this is on inv, all others should be?
        std::cout << "edge in: " << edge_in.edge_id << "translated id: " << edge_in.translated_edge_id << " path to edge in: " << std::endl;
        for (auto e:edge_in.path_from_center){
            result.push_back(e);
            std::cout << e << std::endl;
        }
        std::cout << "edge out: " << edge_out.edge_id << "translated id: " << edge_out.translated_edge_id << " path to edge out: " << std::endl;
        for (auto e:edge_out.path_from_center){
            result.push_back(e);
            std::cout << e << std::endl;
        }
        result.push_back(edge_out.translated_edge_id);
    } else { // really not sure this is correct
        result.push_back(edge_in.translated_edge_id);
        for (auto e:edge_in.path_from_center){
            result.push_back(involution[e]);
        }
        for (auto e:edge_out.path_from_center){
            result.push_back(involution[e]);
        }
        result.push_back(edge_out.translated_edge_id);
    }

    return result;

}


ComplexRegionCollection::ComplexRegionCollection(vec<int>& involution): involution(involution){}



bool ComplexRegionCollection::AddRegion(std::vector<uint64_t> edges_in, std::vector<uint64_t> edges_out,
               vec<int> &involution, int insert_size = 5000){
    std::set<uint64_t > check_edges_distinct;
    for (auto edge: edges_in){
        //std::cout << "adding edge: " << edge << std::endl;
        check_edges_distinct.insert(edge);
        check_edges_distinct.insert(involution[edge]);
        //std::cout << "size of set: " << check_edges_distinct.size() << std::endl;
    }
    for (auto edge: edges_out){
        //std::cout << "adding edge: " << edge << "size of set: " << check_edges_distinct.size() << std::endl;
        check_edges_distinct.insert(edge);
        check_edges_distinct.insert(involution[edge]);
        //std::cout << "size of set: " << check_edges_distinct.size() << std::endl;

    }
    if (check_edges_distinct.size() == (2*edges_in.size() + 2*edges_out.size())) {
        //std::cout << "creating region " << std::endl;
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
        std::cout << "Considering region "<< path_str(region.edges_in) << std::endl;
        if (region.solved){
            solved_regions.push_back(region);
        }
    }
    for (int i = 0; i < solved_regions.size() ; i++){
        auto region = solved_regions[i];
        bool region_clashed = false;
        auto edges_in = region.edges_in;
        std::cout << "edges in: " << path_str(edges_in) << " i: " << i << std::endl;
        for (int j = i + 1; j < solved_regions.size(); j++){
            auto region_2 = solved_regions[j];
            for (auto edge: edges_in) {
                std::cout << "Looking for edge: " << edge << " j: " << j << std::endl;
                if (std::find(region_2.edges_in.begin(), region_2.edges_in.end(), edge) != region_2.edges_in.end()
                        || std::find(region_2.edges_in.begin(), region_2.edges_in.end(), involution[edge]) != region_2.edges_in.end()){
                        std::cout << "edge found: " << path_str(region_2.edges_in) << std::endl;
                        solved_regions_final.push_back(FindBestSolvedRegion(region, region_2));
                        region_clashed = true;
                    std::cout << "region chosen:" << path_str(solved_regions_final[solved_regions_final.size()- 1].edges_in) << std::endl;
                }
            }
        }
        if (!region_clashed){
            //even through this gets hit by a breakpoint, solved regions final size doesn't go up
            solved_regions_final.push_back(region);
            std::cout << "Added solved rgeion: " << solved_regions_final.size() << std::endl;
        }
    }
    std::cout << "Added: " << solved_regions_final.size() << "solved regions"  <<std::endl;
    if (solved_regions_final.size() > 1) {
        for (int i = 0; i < solved_regions_final.size(); i++) {
            bool region_removed = false;
            auto region = solved_regions_final[i];
            auto edges_out = region.edges_out;
            std::cout << "checking edges out: " << path_str(edges_out) << " number solved regions:" <<  solved_regions_final.size() << " i: " << i << " i+1: " << i+1 << std::endl;
            for (int j = i + 1; j < solved_regions_final.size() -1; j++) { // ok somehow j gets to 262149
                auto region_2 = solved_regions_final[j];
                for (auto edge: edges_out) {
                    if (std::find(region_2.edges_out.begin(), region_2.edges_out.end(), edge) != region_2.edges_out.end()
                        || std::find(region_2.edges_out.begin(), region_2.edges_out.end(), involution[edge]) !=
                           region_2.edges_out.end()) {
                        std::cout << "edge found: " << path_str(region_2.edges_in) << std::endl;
                        if (region.edges_in.size() > region_2.edges_in.size()){ //TODO define proper = operator so we cna ue the select best region function here tooo
                            std::cout << "Removing: " << path_str(solved_regions_final[j].edges_out) << std::endl;
                            solved_regions_final.erase(solved_regions_final.begin() + j);
                        } else {
                            std::cout << "Removing: " << path_str(solved_regions_final[i].edges_out) << std::endl;
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
    for (auto region:solved_regions_final){
        std::cout << "Solved region edges in: " << path_str(region.edges_in) << std::endl;
        std::cout << "Solved region edges out: " << path_str(region.edges_out) << std::endl;
    }
    std::cout << "Number of potential solved regions after removing clashing ones: " << solved_regions_final.size() << std::endl;
}


std::vector<std::vector<uint64_t> > ComplexRegionCollection::GetPathsToSeparate(){
    // for each solved region, get each path in the right direction (so from in to out), on the main graph, not the involution
    std::vector<std::vector<uint64_t> > paths_to_separate;
    for (auto region:solved_regions_final){
        std::cout << "trying to solve region" << std::endl; // this is printed out, so issue is below here-
        auto paths = region.BuildPaths();
        for (auto path:paths){
            paths_to_separate.push_back(path);
        }
    }
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

// instead of doing it like this, find pairs of clashing regions and select one which solves most paths
ComplexRegion ComplexRegionCollection::FindBestSolvedRegion(ComplexRegion region_1, ComplexRegion region_2){
    if (region_1.edges_in.size() > region_2.edges_in.size()){
        return region_1;
    } else {
        return region_2;
    }
}
