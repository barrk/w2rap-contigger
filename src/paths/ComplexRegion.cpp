//
// Created by Katie Barr (EI) on 15/12/2016.
//

#include "ComplexRegion.h"


ComplexRegion::ComplexRegion(){};

ComplexRegion::ComplexRegion(std::vector<uint64_t  > edges_in, std::vector<uint64_t  > edges_out, vec<int>& involution, int insert_size=5000):
    edges_in(edges_in), edges_out(edges_out), insert_size(insert_size), involution(involution)
{
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

// actually, this would be too late, these can be checked earlier- may still ned to do other things to validate after path is built
bool ComplexRegion::SanityCheckPath(std::vector<uint64_t> path){
    // sanity check path, ensure its not a cycle (regions with equal in/out edges shouldn't be added), and is in this region
    if (path.front() != path.back() &&
            std::find(edges_in_canonical.begin(), edges_in_canonical.end(), path.front()) != edges_in_canonical.end()
            && std::find(edges_out_canonical.begin(), edges_out_canonical.end(), path.back()) != edges_out_canonical.end()){
        return true; // the path is not a cycle, it starts and ends at the bounds of this region- ask Bernardo what else I should check
    } else {
        std::cout << "Path is not valid in this region" << std::endl;
        return false;
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
            std::cout << "edge:" << edge << " is canonical, involution: " << involution[edge] << std::endl;
            edge_details.edge_id = edge;
            edge_details.forward = true;
            edges_canonical[i] = edge;
        } else {
            std::cout << "edge:" << edge << " canonical version " << edge<< std::endl;
            edge_details.edge_id = involution[edge];
            edge_details.forward = false;
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
            in_edges_solved.insert(count.first.first);
            out_edges_solved.insert(count.first.second);
            number_solved_pairs += 1;
            std::cout << "solved" << std::endl;
        }
    }
    if(in_edges_solved.size() == edges_in.size() && out_edges_solved.size() == edges_out.size()){
        solved = true;
        std::cout << "Region solved" << std::endl;
    }
}

std::vector<std::vector<uint64_t> >  ComplexRegion::BuildPaths(){
    /*
     * determine which in/out pairs to use and build th associated paths
    */
    std::vector<std::vector<uint64_t> > paths;
    // combination counts contains pairs of in/out edge ids,
    for (auto in_out_pair: combination_counts){
        auto in_edge_id = in_out_pair.first.first;
        auto out_edge_id = in_out_pair.first.second;
        std::cout << "Looing for path between " << in_edge_id << " and " << out_edge_id << std::endl;
        for (auto in_edge: edges_in_detailed){
            if (in_edge.edge_id == in_edge_id) {
                for (auto out_edge: edges_in_detailed) {
                    if (out_edge.edge_id == out_edge_id){
                        auto path = BuildPath(in_edge, out_edge);
                        paths.push_back(path);
                    }
                }
            }
        }

    }
    return paths;
}
// this is baed on paths to original edges in/out
std::vector<uint64_t>  ComplexRegion::BuildPath(BoundingEdge edge_in, BoundingEdge edge_out){
    std::cout << "builng path between: "<< edge_in.edge_id <<" and " <<edge_out.edge_id << std::endl;
    std::vector<uint64_t> result;
    if (edge_in.forward){
        result.push_back(edge_in.translated_edge_id); // should be same as edge id- if this is on inv, all others should be?
        for (auto e:edge_in.path_from_center){
            result.push_back(e);
        }
        for (auto e:edge_out.path_from_center){
            result.push_back(e);
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



void ComplexRegionCollection::AddRegion(std::vector<uint64_t> edges_in, std::vector<uint64_t> edges_out,
               vec<int> &involution, int insert_size = 5000){
    std::set<uint64_t > check_edges_distinct;
    for (auto edge: edges_in){
        check_edges_distinct.insert(edge);
    }
    for (auto edge: edges_out){
        check_edges_distinct.insert(edge);
    }
    if (check_edges_distinct.size() == (edges_in.size() + edges_out.size())) {
        ComplexRegion complex_region(edges_in, edges_out, involution, insert_size);
        complex_regions.push_back(complex_region);
        auto key = std::make_pair(complex_region.edges_in, complex_region.edges_out);
        edges_to_region_index[key] = complex_regions.size();
    } else {
        std::cout << "In/Out edges do not define a valid region" << std::endl;
    }
}

bool ComplexRegionCollection::ContainsRegionWithEdges(std::vector<uint64_t> edges_in, std::vector<uint64_t> edges_out){
    return edges_to_region_index.count(std::make_pair(edges_in, edges_out)) != 0;
}

ComplexRegion ComplexRegionCollection::GetRegionWithEdges(std::vector<uint64_t> edges_in, std::vector<uint64_t> edges_out){
    auto index = edges_to_region_index[std::make_pair(edges_in, edges_out)];
    return complex_regions[index];

}


void ComplexRegionCollection::SelectRegionsForPathSeparation(){
    std::vector<std::vector<uint64_t > > solved_region_in;
    std::vector<std::vector<uint64_t > > solved_region_out;
    for (auto region:complex_regions){
        if (region.solved){
            solved_regions.push_back(region);
            solved_region_in.push_back(region.edges_in_canonical);
            solved_region_out.push_back(region.edges_out_canonical);
        }
    }
    std::cout << "Number of potential solved regions: " << solved_regions.size() << std::endl;
    auto distinct_in_edge_sets = CheckNoPathsClash(solved_region_in);
    auto distinct_out_edge_sets = CheckNoPathsClash(solved_region_out);
    if (distinct_in_edge_sets != solved_region_in.size() || distinct_out_edge_sets != solved_region_out.size()){
        std::cout << "Region clash!!" << std::endl; // i really hope this doesn't happen
        // just remove the one that clashes, until none clash
        while (distinct_in_edge_sets != solved_region_in.size()){
            solved_region_in.erase(solved_region_in.begin() + distinct_in_edge_sets);
            solved_region_out.erase(solved_region_out.begin() + distinct_in_edge_sets);
            solved_regions.erase(solved_regions.begin() + distinct_in_edge_sets);
            distinct_in_edge_sets = CheckNoPathsClash(solved_region_in);
        }
        while (distinct_out_edge_sets != solved_region_out.size()){
            solved_region_in.erase(solved_region_in.begin() + distinct_out_edge_sets);
            solved_region_out.erase(solved_region_out.begin() + distinct_out_edge_sets);
            solved_regions.erase(solved_regions.begin() + distinct_out_edge_sets);
            distinct_out_edge_sets = CheckNoPathsClash(solved_region_out);
        }
    }
    std::cout << "Number of potential solved regions after removing clashing ones: " << solved_regions.size() << std::endl;
}


std::vector<std::vector<uint64_t> > ComplexRegionCollection::GetPathsToSeparate(){
    // for each solved region, get each path in the right direction (so from in to out), on the main graph, not the involution
    std::vector<std::vector<uint64_t> > paths_to_separate;
    for (auto region:solved_regions){
        auto paths = region.BuildPaths();
        for (auto path:paths){
            paths_to_separate.push_back(path);
        }
    }
    return paths_to_separate;
}



int ComplexRegionCollection::CheckNoPathsClash(std::vector<std::vector<uint64_t > > all_edges){
    std::set<uint64_t > seen_edges;

    auto seen_count = seen_edges.size();
    int region_index = 0;
    for (auto edges: all_edges){
        for (auto edge:edges){
            seen_edges.insert(edge);
            if (seen_edges.size() == (seen_count + 1)){
                seen_count += 1;
            } else{
                return region_index;
            }
        }
        region_index += 1;
    }
    // this gives us the index at which one solved region clashes with another
    return region_index;
}