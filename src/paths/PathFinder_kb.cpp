//
// Created by Katie Barr (EI) on 31/10/2016.
//
//


#include "PathFinder_kb.h"
using namespace PatherKB;

void PathFinderkb::init_prev_next_vectors(){
    //TODO: this is stupid duplication of the digraph class, but it's so weird!!!
    prev_edges.resize(mToLeft.size());
    next_edges.resize(mToRight.size());
    for (auto e=0;e<mToLeft.size();++e){

        uint64_t prev_node=mToLeft[e];

        prev_edges[e].resize(mHBV.ToSize(prev_node));
        for (int i=0;i<mHBV.ToSize(prev_node);++i){
            prev_edges[e][i]=mHBV.EdgeObjectIndexByIndexTo(prev_node,i);
        }

        uint64_t next_node=mToRight[e];

        next_edges[e].resize(mHBV.FromSize(next_node));
        for (int i=0;i<mHBV.FromSize(next_node);++i){
            next_edges[e][i]=mHBV.EdgeObjectIndexByIndexFrom(next_node,i);
        }
    }


}


std::tuple <std::vector<LMPPair >, std::vector<LMPPair >, std::map<uint64_t, std::vector<int> > > PathFinderkb::mapEdgesToLMPReads(){
    KMatch kmatch(31);
    LMPMapper lmp_mapper(lmp_data, mHBV, mInv, kmatch);
    std::vector<LMPPair >  lmp_pairs;
    std::map<uint64_t, std::vector<int> > edge_id_to_pair_id_map;
    std::vector<LMPPair >  lmp_pairs_for_insert_size_estimation;
    lmp_mapper.LMPReads2MappedPairedEdgePaths(lmp_pairs, lmp_pairs_for_insert_size_estimation, edge_id_to_pair_id_map);
    std::cout << "mapEdgesToLMPReads lmp pairs size: " << lmp_pairs.size() << std::endl;
    return std::make_tuple(lmp_pairs, lmp_pairs_for_insert_size_estimation, edge_id_to_pair_id_map);
}


// this allows it to go round in circles indefinitely, which we don't want
void PathFinderkb::edges_beyond_distance(std::vector<uint64_t>  & spanning_edges, std::vector<std::vector<uint64_t> >  & paths_to_spanning_edges, std::vector<uint64_t> & intermediate_path, uint64_t e, std::vector<uint64_t > & traversed_edge_list, uint64_t approximate_insert_size, int recursion_depth=0, int distance_traversed=0, std::string direction="right") {
    /*
     * to find edges in region lmp pairs might solve, we want to traverse graph until we're large_frontier_size away from e
        looking at some of the thins returned from this for larger test datasets, i'm not covinced it does exactly what we want
        as it quite often returns lots and lots of paths
        maybe if we did something like detected if we've traversed several forks, which would make the region too coplex to solve,
        or moved ont othe RC graph, we could have it return something to indicate that this region is too complex
        though to be fair, returning a massive list does indicate that!!
        for time being just return edges which are long enough to be exiting a complex region- like bernardo's long frontier size
     */
    std::vector<uint64_t> edges;
    if (direction=="right"){edges = prev_edges[e];

    } else {
        edges = next_edges[e];

    }

    // edges between the start edge and the end edge are ones we traverse through in that recursion stack,
    // only want paths from beginning node, i.e. when recursion depth is 0, but paths can split at any time
    for (auto edge:edges){
        auto edge_length = mHBV.EdgeObject(edge).size();
        // if we haven't traversed this edge, traverse it
        if (std::find(traversed_edge_list.begin(), traversed_edge_list.end(), edge) == traversed_edge_list.end()) {
            // if this edge takes us far enough away, add it to long fronteirs
            if ((edge_length > (approximate_insert_size - distance_traversed))){// && (edge_length > 10000)){ segfaults if i add this!
                spanning_edges.push_back(edge);
                // i want to store edges between start edge and edge, this stores attached edges, stores the ones traversed before twice for some readon
                traversed_edge_list.push_back(edge);
                std::vector<uint64_t> temp_path;

                for (auto intermediate_edge: intermediate_path){
                    temp_path.push_back(intermediate_edge);
                }
                paths_to_spanning_edges.push_back(temp_path);
                temp_path.clear();
            } else { // if this edge does not take us far enough away, add its length to distance traversed go to next edge
                if (direction=="left"){
                    intermediate_path.insert(intermediate_path.begin(), edge);
                } else {
                    intermediate_path.push_back(edge);
                }
                distance_traversed += mHBV.EdgeObject(edge).size();
                recursion_depth += 1;
                edges_beyond_distance(spanning_edges, paths_to_spanning_edges, intermediate_path, edge, traversed_edge_list, approximate_insert_size,  recursion_depth, distance_traversed, direction);
                intermediate_path.clear();
            }
        }
    }

}

std::vector<uint64_t>  PathFinderkb::canonicalisePath(std::vector<uint64_t> path){
    std::vector<uint64_t> result;
    bool involution_path = false;
    // if first edge is on involution, all will be
    if (path[0] > mInv[path[0]]){
        involution_path = true;
    }
    for (auto edge: path){
        if (involution_path){
            result.push_back(mInv[edge]);
        } else {
            result.push_back(edge);
        }
    }
    //if (result.back() < result.front()){
    //    std::reverse(result.begin(), result.end());
    //}
    if (involution_path || result.back() < result.front()) {
        std::reverse(result.begin(), result.end());
    }
    return result;

}



typedef struct {
    uint64_t e1;
    uint64_t e2;
    int path_id;
    int length;
} PathDetails;

void PathFinderkb::resolveComplexRegionsUsingLMPData() {

    int approximate_insert_size = 3500;
    init_prev_next_vectors();// mToLeft and mToright are to/from vertices
    std::tuple<std::vector<LMPPair>, std::vector<LMPPair>, std::map<uint64_t, std::vector<int> > > mapping_results = mapEdgesToLMPReads();
    std::vector<LMPPair> pairs_for_scaffolding = std::get<0>(mapping_results);
    std::vector<LMPPair> pairs_for_insert_size_calculation = std::get<1>(mapping_results);
    std::map<uint64_t, std::vector<int> > edge_id_to_pair_id_map = std::get<2>(mapping_results);
    int same_in_out_degree = 0;
    int same_in_out_degree_complex = 0;
    int solveable_regions_count = 0;

    vector<uint64_t> traversed_edge_list;
    std::vector<std::vector<uint64_t> > paths_to_spanning_edges;
    std::set<std::vector<uint64_t> > paths_seen;
    std::vector<uint64_t> intermediate_path;
    std::vector<uint64_t> involutions_of_edges_in_solved_regions;
    std::vector<uint64_t> spanning_edges_in;
    std::vector<uint64_t> spanning_edges_out;
    ComplexRegion complex_region;
    ComplexRegionCollection complex_regions(mInv);

    //std::vector<int> edges = {8, 9, 320, 390, 336, 463};
    for (int edge_index = 0; edge_index < mHBV.EdgeObjectCount(); ++edge_index) {
        // if we have already seen the involution of this edge, then leave it
        if (std::find(involutions_of_edges_in_solved_regions.begin(), involutions_of_edges_in_solved_regions.end(), edge_index) == involutions_of_edges_in_solved_regions.end()) {
            auto edge = mHBV.EdgeObject(edge_index).ToString();
            if (edge == mHBV.EdgeObject(mInv[edge_index]).ToString()) {
                continue;
            }
            edges_beyond_distance(spanning_edges_in, paths_to_spanning_edges, intermediate_path, edge_index,
                                  traversed_edge_list, approximate_insert_size, 0, 0, "right");
            traversed_edge_list.clear();
            edges_beyond_distance(spanning_edges_out, paths_to_spanning_edges, intermediate_path, edge_index,
                                  traversed_edge_list, approximate_insert_size, 0, 0, "left");
            //std::cout << "in" << path_str(spanning_edges_in) << std::endl;
            //std::cout << "out " <<  path_str(spanning_edges_out) << std::endl;
            // take regions with less than 5 in/out edges- ideally would make graph traversal give up if it finds too many
            if ((spanning_edges_in.size() > 1) && (spanning_edges_in.size() == spanning_edges_out.size()) &&
                spanning_edges_in.size() < 5) {
                std::cout << "in" << path_str(spanning_edges_in) << std::endl;
                std::cout << "out " << path_str(spanning_edges_out) << std::endl;
                std::cout << "complex regions size: " << complex_regions.complex_regions.size() << std::endl;
                if (complex_regions.ContainsRegionWithEdges(spanning_edges_in, spanning_edges_out)) {
                    // don't think this should happen, but we may have overlapping regions where we just want to selet one
                    complex_region = complex_regions.GetRegionWithEdges(spanning_edges_in, spanning_edges_out);
                } else {// if we have't already created this region, do so
                    // try to move this further down- or keep index from beginning- i.e. add counter to for
                    bool region_added = complex_regions.AddRegion(spanning_edges_in, spanning_edges_out, mInv,
                                                                  approximate_insert_size);
                    if (region_added) {
                        std::cout << "comple regions size" << complex_regions.complex_regions.size() << std::endl;
                        complex_region = complex_regions.complex_regions.back();
                        std::cout << "egdes in of region added: " << path_str(complex_region.edges_in) << std::endl;
                        std::cout << "edges out of region added: " << path_str(complex_region.edges_out) << std::endl;
                        std::cout << "comple regions size" << complex_regions.complex_regions.size() << std::endl;
                    } else {
                        continue;
                    }
                }
                std::cout << "complex region created" << std::endl;
                same_in_out_degree_complex += 1;
                int count = 0;
                int edge_ind = 0;
                for (auto edge_in: spanning_edges_in) {
                    std::vector<uint64_t> intermediate_edges;
                    for (auto int_edge: paths_to_spanning_edges[count]) {
                        intermediate_edges.push_back(int_edge);
                    }
                    intermediate_edges.push_back(edge_index);
                    //std::cout << "Path for edge: " << complex_region.edges_in[edge_ind] << " " <<  path_str(intermediate_edges) << std::endl;
                    complex_region.AddPathTo(edge_ind, true, intermediate_edges);
                    //std::cout << "Ids of mapping pairs: " << path_str(edge_id_to_pair_id_map[edge_in]) << std::endl;
                    for (auto pair_id :  edge_id_to_pair_id_map[edge_in]) {
                        complex_region.AddPairId(edge_ind, pair_id, true);
                    }
                    count += 1;
                    edge_ind += 1;
                }
                edge_ind = 0;
                for (auto edge_out: spanning_edges_out) {
                    std::vector<uint64_t> intermediate_edges;
                    for (auto int_edge: paths_to_spanning_edges[count]) {
                        intermediate_edges.push_back(int_edge);
                    }
                    complex_region.AddPathTo(edge_ind, false, intermediate_edges);
                    //std::cout << "Path for edge: " << complex_region.edges_in[edge_ind] << " " <<  path_str(intermediate_edges) << std::endl;
                    //std::cout << "Ids of mapping pairs: " << path_str(edge_id_to_pair_id_map[edge_out]) << std::endl;
                    for (auto pair_id :  edge_id_to_pair_id_map[mInv[edge_out]]) {
                        complex_region.AddPairId(edge_ind, pair_id, false);
                    }
                    count += 1;
                    edge_ind += 1;

                }
                complex_region.FindSolveablePairs();
                complex_region.isSolved(5);
                if (complex_region.solved) {
                    std::cout << "egdes in of region solved: " << path_str(complex_region.edges_in) << std::endl;
                    std::cout << "edges out of region solved: " << path_str(complex_region.edges_out) << std::endl;
                    std::cout << "complex regions size: " << complex_regions.complex_regions.size() << std::endl;
                    solveable_regions_count += 1;
                    involutions_of_edges_in_solved_regions.push_back(mInv[edge_index]);
                }
                complex_regions.complex_regions.back() = complex_region;
            }
            if (edge_index ==321 || edge_index ==322 || edge_index ==323 || edge_index ==324 || (edge_index > mHBV.EdgeObjectCount() -3)) {

                for (auto region: complex_regions.complex_regions) {
                    std::cout << "during loop, edge index " << edge_index <<" region with in edges: " << path_str(region.edges_in)
                              << " solved boolean " << region.solved << std::endl;
                }
            }
            paths_to_spanning_edges.clear();
            spanning_edges_in.clear();
            spanning_edges_out.clear();
            /*
            if (mInv[edge_index + 1] == edge_index){
                edge_index += 1; // if the next edge will be involution of current one, skip it
            } = actually better to properly check */
        }
    }
    for (auto region: complex_regions.complex_regions){
        std::cout << "region with in edges: " << path_str(region.edges_in) << " solved boolean " <<region.solved << std::endl;
    }
    std::cout << "Solveable regions: " << solveable_regions_count << " of: " << complex_regions.complex_regions.size() << std::endl;
    // can now compare regions, select ones which are solved, track paths to ensure ends don't meet
    complex_regions.SelectRegionsForPathSeparation();
    auto paths_to_separate = complex_regions.GetPathsToSeparate();
    // last bit hsould be same as before
    uint64_t sep=0;
    std::map<uint64_t,std::vector<uint64_t>> old_edges_to_new;
    for (auto p:paths_to_separate){
        std::cout << "separating path: " << path_str(p) << std::endl;
        if (old_edges_to_new.count(p.front()) > 0 or old_edges_to_new.count(p.back()) > 0) {
            std::cout<<"WARNING: path starts or ends in an already modified edge, skipping"<<std::endl;
            continue;
        }
            auto oen=separate_path(p, true);
            if (oen.size() > 0) {
                for (auto et:oen) {
                    //std::cout << "Edge to migrate: " << et.first;
                        if (old_edges_to_new.count(et.first) == 0) old_edges_to_new[et.first] = {};
                        for (auto ne:et.second) {
                            old_edges_to_new[et.first].push_back(ne);
                            //std::cout << " " << ne;
                        }
                    //std::cout << std::endl;
                }
                sep++;
            }
        }
    //mHBV.Involution(mInv);
    //TestInvolution(mHBV, mInv);

    // pretty sure, because these are lmp paths,
    //if (old_edges_to_new.size()>0) {
    //    migrate_readpaths(old_edges_to_new);
    //}
    std::cout<<" "<<sep<<" paths separated!"<<std::endl;
    std::cout<<mHBV.EdgeObjectCount() << " edges in hbv"<<std::endl;
    std::cout << "Inv edge 0: " << mInv[0] << " edge 494: " << mInv[494] << " edge 495 " << mInv[495] << " edge 1: " << mInv[1] << " edge 8 " << mInv[8] << " edge 9: " << mInv[9] << std::endl;
    TestInvolution(mHBV, mInv);
    mHBV.Involution(mInv);
    std::cout << "Inv edge 0: " << mInv[0] << " edge 494: " << mInv[494] << " edge 495 " << mInv[495] << " edge 1: " << mInv[1] << " edge 8 " << mInv[8] << " edge 9: " << mInv[9] << std::endl;

    TestInvolution(mHBV, mInv);
}

void PathFinderkb::resolveRegionsUsingLMPData() {

    int approximate_insert_size = 3500;
    init_prev_next_vectors();// mToLeft and mToright are to/from vertices
    std::tuple<std::vector<LMPPair >, std::vector<LMPPair >, std::map<uint64_t, std::vector<int> > > mapping_results = mapEdgesToLMPReads();
    std::vector<LMPPair > pairs_for_scaffolding = std::get<0>(mapping_results);
    std::vector<LMPPair > pairs_for_insert_size_calculation = std::get<1>(mapping_results);
    std::map<uint64_t, std::vector<int> >  edge_id_to_pair_id_map = std::get<2>(mapping_results);
    int same_in_out_degree = 0;
    int same_in_out_degree_complex = 0;
    int  solveable_regions_count = 0 ;
    std::map<uint64_t, int> mapping_counts;

    vector<uint64_t > traversed_edge_list;
    std::vector<std::vector<uint64_t> >  paths_to_spanning_edges;
    std::set<std::vector<uint64_t> >  paths_seen;
    std::vector<uint64_t> intermediate_path;
    std::vector<uint64_t>  spanning_edges_in;
    std::vector<uint64_t>  spanning_edges_out;
    std::vector<std::vector<uint64_t> > paths_to_separate;

    std::vector<int> edges = {320, 390, 336, 463};
    //for (int edge_index = 0; edge_index < mHBV.EdgeObjectCount(); ++edge_index) {
    for (int edge_index:edges) {
    //for (int edge_index = 0; edge_index < mHBV.EdgeObjectCount(); ++edge_index) {
        auto edge = mHBV.EdgeObject(edge_index).ToString();
        if (edge == mHBV.EdgeObject(mInv[edge_index]).ToString()){
            continue;
        }
            edges_beyond_distance(spanning_edges_in, paths_to_spanning_edges, intermediate_path, edge_index, traversed_edge_list, approximate_insert_size, 0, 0, "right");
            traversed_edge_list.clear();
            edges_beyond_distance(spanning_edges_out, paths_to_spanning_edges, intermediate_path, edge_index, traversed_edge_list, approximate_insert_size, 0, 0, "left");

            if (spanning_edges_in.size() == spanning_edges_out.size()) {//} && long_frontiers_in.size() != 0) {
                same_in_out_degree += 1;
                if (spanning_edges_in.size() > 1) {
                    same_in_out_degree_complex += 1;
                    std::map<uint64_t, std::pair< std::vector<uint64_t>, std::vector<int> > > mapped_lmp_in;
                    std::map<uint64_t, std::pair< std::vector<uint64_t>, std::vector<int> > > mapped_lmp_out;
                    int count = 0;
                    // find number of lmp reads mapping to each large fronteir edge
                    for (auto edge_in: spanning_edges_in) {
                        std::vector<uint64_t> intermediate_edges;
                        for (auto int_edge: paths_to_spanning_edges[count]){
                            intermediate_edges.push_back(int_edge);
                        }
                        std::vector<int> pair_ids;
                        for (auto pair_id :  edge_id_to_pair_id_map[edge_in]) {
                            pair_ids.push_back(pair_id);
                            mapping_counts[edge_in] += 1;
                        }
                        count += 1;
                        std::pair< std::vector<uint64_t>, std::vector<int> > result = std::make_pair(intermediate_edges, pair_ids);
                        mapped_lmp_in[edge_in] = result;
                    }
                    for (auto edge_out: spanning_edges_out) {
                        std::vector<uint64_t> intermediate_edges;
                        for (auto int_edge: paths_to_spanning_edges[count]){
                            intermediate_edges.push_back(int_edge);
                        }
                        std::vector<int> pair_ids;
                        for (auto pair_id :  edge_id_to_pair_id_map[mInv[edge_out]]) {
                            mapping_counts[edge_out] += 1;
                            pair_ids.push_back(pair_id);

                        }
                        count += 1;
                        std::pair< std::vector<uint64_t>, std::vector<int> > result = std::make_pair(intermediate_edges, pair_ids);
                        mapped_lmp_out[edge_out] = result;

                    }
                    std::map<std::pair<uint64_t, uint64_t>, int > counts_for_each_solveabe_pair;
                    // find pairs which are in both mapped_lmp_in, and mapped_lmp_out- these can be used to solve this region
                    for (auto edge_in: spanning_edges_in) {
                        // TODO this way is horrible, need to design a better data structure
                        vector<int> p1_ids = mapped_lmp_in[edge_in].second;
                        for (auto p1_id: p1_ids){
                            for (auto edge_out: spanning_edges_out) {
                                vector<int> p2_ids = mapped_lmp_out[edge_out].second;
                                if ((std::find(p2_ids.begin(), p2_ids.end(), p1_id) !=
                                        p2_ids.end())) {
                                    counts_for_each_solveabe_pair[std::make_pair(edge_in, edge_out)] += 1;

                                }
                            }
                        }
                    }

                    int paths_separated = 0;
                    std::vector<PathDetails> path_end_tracker;
                    for (auto edge_pair_count: counts_for_each_solveabe_pair){
                        std::pair<uint64_t, uint64_t> key = edge_pair_count.first;
                        uint64_t edge_in = key.first;
                        uint64_t edge_out = key.second;
                        int count = edge_pair_count.second;
                        std::cout << "count for pair:" << edge_pair_count.first.first <<  " " << edge_pair_count.first.second << " : " <<count <<std::endl;
                        if (count > 10){
                            std::vector<uint64_t> path_in = mapped_lmp_in[edge_in].first;
                            std::vector<uint64_t> path_out = mapped_lmp_out[edge_out].first;
                            std::vector<uint64_t> full_path;
                            full_path.push_back(edge_in);
                            //std::cout << "edge in added;" << path_str(full_path) << std::endl;
                            for (auto e:path_in){
                                full_path.push_back(e);
                            }
                            //std::cout << "path in added;" << path_str(full_path) << std::endl;
                            full_path.push_back(edge_index);
                            //std::cout << path_str(full_path) << std::endl;
                            for (auto e:path_out){
                                full_path.push_back(e);
                            }
                            //std::cout << "edge;" << path_str(full_path) << std::endl;
                            full_path.push_back(edge_out);
                            std::cout << "path out added;" << path_str(full_path) << std::endl;
                            std::vector<uint64_t> path_canonical = canonicalisePath(full_path);
                            std::vector<uint64_t> path_reversed;
                            path_reversed.resize(path_canonical.size());
                            std::reverse_copy(path_canonical.begin(), path_canonical.end(), path_reversed.begin());
                            std::vector<uint64_t> edges_to_look_for {77, 76, 197, 196};
                            if ((std::find(paths_seen.begin(), paths_seen.end(), path_canonical) == paths_seen.end()) && (std::find(paths_seen.begin(), paths_seen.end(), path_reversed) == paths_seen.end())) {
                                int inner_path_length = 0;
                                for (int i =1; i < path_canonical.size() -1; i++){
                                    auto edge = mHBV.EdgeObject(path_canonical[i]);
                                    inner_path_length += edge.size();
                                }
                                if (inner_path_length < 12000) { // hardcoded something longer tha any lmp... insert size
                                    paths_seen.insert(path_canonical);
                                    std::cout << "path to separate: " << path_str(full_path) << std::endl;
                                    std::cout << "path canonical: " << path_str(path_canonical) << std::endl;
                                    paths_to_separate.push_back(path_canonical);
                                    paths_separated += 1;
                                }
                            }
                        }
                    }
                    std::cout << "total paths eparated: " << paths_separated << " out of " << spanning_edges_in.size() << " len paths to separate " << paths_to_separate.size() << std::endl;
                    if (spanning_edges_in.size() == paths_separated){
                        solveable_regions_count += 1;
                    }
                }

            }
            paths_to_spanning_edges.clear();
            spanning_edges_in.clear();
            spanning_edges_out.clear();
    }
    // get both [354:355 204:205 481:480 ] and [354:355 204:205 302:303 204:205 481:480 ], want the longer- these must occur for different edges, so need tobe out of main loop
    // should be able to incorporate this into above somewhere- but need to know length of path. given the paths are already calculated that can probably be moved up too
    /*std::set<std::set<uint64_t, uint64_t> > edge_pairs;
    for (std::vector<uint64_t> path; paths_to_separate){
        auto before = edge_pairs.size();
        std::set<uint64_t, uint64_t> path_ends = {path.front(), path.back()};
        edge_pairs.insert(path_ends);
        auto after = edge_pairs.size();
        if (before == after){
            // this doesn't do it, because we still need to know length of path, and index of path with these ends

        }


    }*/

    std::map<std::pair<uint64_t, uint64_t>, PathDetails> path_length_edge_map;
    int path_index = 0;
    for (auto path: paths_to_separate){
        auto path_length = path.size();
        PathDetails path_details;
        path_details.path_id = path_index;
        path_details.length = path_length;
        uint64_t start = path.front();
        uint64_t end = path.back();
        uint64_t start_rc = mInv[start];
        uint64_t end_rc = mInv[end];
        std::pair<uint64_t, uint64_t> edges;
        edges.first = start;
        edges.second = end;
        std::pair<uint64_t, uint64_t>  edges_rc;
        edges_rc.first = start_rc;
        edges_rc.second = end_rc;
        //std::cout << "edge in:" << start << " edge out:" << end <<std::endl;
        std::cout << path_length_edge_map.count(edges) << std::endl;
        //std::cout << "edge in rc:" << start_rc << " edge out rc:" << end_rc <<std::endl;
        std::cout << path_length_edge_map.count(edges_rc) << std::endl;
        std::pair<uint64_t, uint64_t> edges_backwards;
        edges_backwards.first = edges.second;
        edges_backwards.second = edges.first;
        std::pair<uint64_t, uint64_t> edges_backwards_rc;
        edges_backwards_rc.first = edges_rc.second;
        edges_backwards_rc.second = edges_rc.first;
        if (path_length_edge_map.count(edges) == 0 && path_length_edge_map.count(edges_rc) == 0 && (path_length_edge_map.count(edges_backwards) == 0 && path_length_edge_map.count(edges_backwards_rc) == 0)){
            path_length_edge_map[edges] = path_details;
            //std::cout << "edges p len: " << path_length_edge_map[edges].length << std::endl;
            //std::cout << "edges_rc p len: " << path_length_edge_map[edges_rc].length << std::endl;
            if ((path_length_edge_map.count(edges) == 0) && (path_length_edge_map.count(edges_rc) == 0)){
                //std::cout <<  " first cond part:" << std::endl;
            }
            //std::cout << path_str(path) << std::endl;
            //std::cout << "count edges: " << path_length_edge_map.count(edges)<< std::endl;
            //std::cout << "count edges_rc: " << path_length_edge_map.count(edges_rc) << std::endl;
            for (auto pmap:path_length_edge_map){
                std::cout << "key: " << pmap.first.first << " " << pmap.first.second << ", value: " << pmap.second.path_id << std::endl;
            }
        }
        path_index += 1;
    }

    std::vector<std::vector<uint64_t> > paths_to_separate_final;
    for (auto path_to_use: path_length_edge_map){
        std::cout << "path str: " << path_to_use.second.path_id << std::endl;
        paths_to_separate_final.push_back(paths_to_separate[path_to_use.second.path_id]);
        std::cout << path_str(paths_to_separate[path_to_use.second.path_id]) << std::endl;

    }

    std::cout << "Regions with same degree in and out:" << same_in_out_degree << std::endl;
    std::cout << "Complex regions with same degree in and out:" << same_in_out_degree_complex << std::endl;
    std::cout << "Complex regions solveable wit lmp reads:" << solveable_regions_count<< std::endl;
    uint64_t sep=0;
    std::map<uint64_t,std::vector<uint64_t>> old_edges_to_new;
    std::set<uint64_t> start_edges_seen;
    std::set<uint64_t> end_edges_seen; // old edges to new should do this!!
    for (auto p:paths_to_separate_final){

        if (old_edges_to_new.count(p.front()) > 0 or old_edges_to_new.count(p.back()) > 0) {
            std::cout<<"WARNING: path starts or ends in an already modified edge, skipping"<<std::endl;
            continue;
        }
        uint64_t start = p.front();
        uint64_t end = p.back();
        uint64_t start_rc = mInv[p.front()];
        uint64_t end_rc = mInv[p.back()];
        // this creates new vertices at the end of the graph, connects start to these, and removes old connection then repeats on the involution
        if (std::find(start_edges_seen.begin(), start_edges_seen.end(), start) == start_edges_seen.end() &&
             std::find(end_edges_seen.begin(), end_edges_seen.end(), end) == end_edges_seen.end() &&
                std::find(start_edges_seen.begin(), start_edges_seen.end(), end_rc) == start_edges_seen.end() &&
                std::find(end_edges_seen.begin(), end_edges_seen.end(), start_rc) == end_edges_seen.end()){
                start_edges_seen.insert(start);
                end_edges_seen.insert(end);
                start_edges_seen.insert(end_rc);
                end_edges_seen.insert(start_rc);
                auto oen=separate_path(p, true);
                    if (oen.size() > 0) {
                        //mHBV.Involution(mInv);
                        //TestInvolution(mHBV, mInv);

                        for (auto et:oen) {
                            if (oen[et.first].size() > 0) { // absutely no idea why the start and end edges are ending up as keys in here, because they aren't there when you output oen in side separate paths
                                if (old_edges_to_new.count(et.first) == 0) old_edges_to_new[et.first] = {};
                                for (auto ne:et.second) {
                                    old_edges_to_new[et.first].push_back(ne);
                                }
                            }
                        }
                        sep++;
            }
        }
    }

    if (old_edges_to_new.size()>0) {
        migrate_readpaths(old_edges_to_new);
    }
    std::cout<<" "<<sep<<" paths separated!"<<std::endl;
    std::cout<<mHBV.EdgeObjectCount() << " edges in hbv"<<std::endl;
    // cleanup needs both PE and LMP paths, but then won't know whether its a PE or LMP read- i don't think this matters
    /*std::vector<std::vector<ReadPath> > lmp_paths;
    for (auto path: paths_to_separate_final){
        std::vector<ReadPath>  path_read_path;
        for (auto edge: path){
            path_read_path.push_back(edge);
        }
        lmp_paths.push_back(path_read_path);

    }
    std::vector<std::vector<uint64_t> > all_paths;
    all_paths.insert(all_paths.begin(), mPaths.begin(), mPaths.end());
    all_paths.insert(all_paths.end(), lmp_paths.begin(), lmp_paths.end());*/
    //mHBV.Involution(mInv);
    //TestInvolution(mHBV, mInv);
    //TODO: work out how to convert these paths into read paths!
    //Cleanup(mHBV, mInv, mPaths);


}

std::string PathFinderkb::path_str(std::vector<uint64_t> path) {
    std::string s="[";
    for (auto p:path){
        // output edge id on hbv and involution
        s+=std::to_string(p)+":"+std::to_string(mInv[p])+" ";//+" ("+std::to_string(mHBV.EdgeObject(p).size())+"bp "+std::to_string(paths_per_kbp(p))+"ppk)  ";
    }
    s+="]";
    return s;
}


std::map<uint64_t,std::vector<uint64_t>> PathFinderkb::separate_path(std::vector<uint64_t> p, bool verbose_separation){

    //TODO XXX: proposed version 1 (never implemented)
    //Creates new edges for the "repeaty" parts of the path (either those shared with other edges or those appearing multiple times in this path).
    //moves paths across to the new reapeat instances as needed
    //changes neighbourhood (i.e. creates new vertices and moves the to and from for the implicated edges).

    //creates a copy of each node but the first and the last, connects  linearly to the previous copy,
    std::cout<<std::endl<<"Separating path"<< path_str(p) << std::endl;
    std::set<uint64_t> edges_fw;
    std::set<uint64_t> edges_rev;
    for (auto e:p){//TODO: this is far too astringent...
        edges_fw.insert(e);
        edges_rev.insert(mInv[e]);

        if (edges_fw.count(mInv[e]) ||edges_rev.count(e) ){ //std::cout<<"PALINDROME edge detected, aborting!!!!"<<std::endl;
            return {};}
    }
    //create two new vertices (for the FW and BW path)
    uint64_t current_vertex_fw=mHBV.N(),current_vertex_rev=mHBV.N()+1; // .N is number of vertices, so these ar ethe two vertices added on the next line
    mHBV.AddVertices(2);
    //migrate connections (dangerous!!!)
    mHBV.GiveEdgeNewToVx(p[0],mToRight[p[0]],current_vertex_fw); // edit graph so edge p now goes to newly created vertex instead of old vertex
    mToRight[p[0]] = current_vertex_fw;
    mHBV.GiveEdgeNewFromVx(mInv[p[0]],mToLeft[mInv[p[0]]],current_vertex_rev);// erases old connection from old involution vertex, and connects edge to old
    mToLeft[mInv[p[0]]] = current_vertex_rev;
    std::map<uint64_t,std::vector<uint64_t>> old_edges_to_new;

    for (auto ei=1;ei<p.size()-1;++ei){
        //add a new vertex for each of FW and BW paths
        uint64_t prev_vertex_fw=current_vertex_fw,prev_vertex_rev=current_vertex_rev;
        //create two new vertices (for the FW and BW path)
        current_vertex_fw=mHBV.N();
        current_vertex_rev=mHBV.N()+1;
        mHBV.AddVertices(2);

        //now, duplicate next edge for the FW and reverse path
        auto nef=mHBV.AddEdge(prev_vertex_fw,current_vertex_fw,mHBV.EdgeObject(p[ei]));// add an edge for each edge in path, with appropriate basevector
        mToLeft.push_back(prev_vertex_fw);
        mToRight.push_back(current_vertex_fw);
        if (! old_edges_to_new.count(p[ei]))  old_edges_to_new[p[ei]]={};
        old_edges_to_new[p[ei]].push_back(nef);

        //std::cout << "Separating: " << p[ei] << " new edge: " << nef << std::endl;
        auto ner=mHBV.AddEdge(current_vertex_rev,prev_vertex_rev,mHBV.EdgeObject(mInv[p[ei]]));
        mToLeft.push_back(current_vertex_rev);
        mToRight.push_back(prev_vertex_rev);
        if (! old_edges_to_new.count(mInv[p[ei]]))  old_edges_to_new[mInv[p[ei]]]={};
        old_edges_to_new[mInv[p[ei]]].push_back(ner);
        //std::cout << "Separating: " << mInv[p[ei]] << " new edge: " << ner << std::endl;

        mInv.push_back(ner);
        mInv.push_back(nef);
        mEdgeToPathIds.resize(mEdgeToPathIds.size()+2);

    }
    mHBV.GiveEdgeNewFromVx(p[p.size()-1],mToLeft[p[p.size()-1]],current_vertex_fw);// attach new edges back into graph at right position
    //mToLeft[p[p.size()-1]] = current_vertex_fw;
    mHBV.GiveEdgeNewToVx(mInv[p[p.size()-1]],mToRight[mInv[p[p.size()-1]]],current_vertex_rev);
    //mToRight[mInv[p[p.size()-1]]] = current_vertex_rev;

    //TODO: cleanup new isolated elements and leading-nowhere paths.
    //for (auto ei=1;ei<p.size()-1;++ei) mHBV.DeleteEdges({p[ei]});
    return old_edges_to_new;

}


void PathFinderkb::migrate_readpaths(std::map<uint64_t,std::vector<uint64_t>> edgemap){
    //Migrate readpaths: this changes the readpaths from old edges to new edges
    //if an old edge has more than one new edge it tries all combinations until it gets the paths to map
    //if more than one combination is valid, this chooses at random among them (could be done better? should the path be duplicated?)
    mHBV.ToLeft(mToLeft);
    mHBV.ToRight(mToRight);
    std::ofstream edge_migrations;
    edge_migrations.open("/Users/barrk/Documents/ecoli_data/v1/edge_migrations.txt");
    for (auto &p:mPaths){
        std::vector<std::vector<uint64_t>> possible_new_edges;
        bool translated=false,ambiguous=false;
        for (auto i=0;i<p.size();++i){
            if (edgemap.count(p[i])) { // if edge i on path p has been translated
                possible_new_edges.push_back(edgemap[p[i]]);
                if (not translated) translated=true;
                // if same edge is in multiple paths, this could occur
                if (possible_new_edges.back().size()>1) ambiguous=true; // if there is more than 1 edge mapping, its ambiguous
            }
            else possible_new_edges.push_back({p[i]});

        }
        if (translated){
            if (not ambiguous){ //just straigh forward translation
                for (auto i=0;i<p.size();++i) {
                    //edge_migrations << "Migrating edge: " << p[i] << " to edge " << possible_new_edges[i][0] << std::endl;
                    p[i]=possible_new_edges[i][0];}
            }
            else {
                //ok, this is the complicated case, we first generate all possible combinations
                std::vector<std::vector<uint64_t>> possible_paths={{}};
                for (auto i=0;i<p.size();++i) {//for each position
                    std::vector<std::vector<uint64_t>> new_possible_paths;
                    for (auto pp:possible_paths) { //take every possible one
                        for (auto e:possible_new_edges[i]) {
                            //if i>0 check there is a real connection to the previous edge
                            if (i == 0 or (mToRight[pp.back()]==mToLeft[e])) {
                                new_possible_paths.push_back(pp);
                                new_possible_paths.back().push_back(e);
                            }
                        }
                    }
                    possible_paths=new_possible_paths;
                    if (possible_paths.size()==0) break;
                }
                if (possible_paths.size()==0){
                    std::cout<<"Warning, a path could not be updated, truncating it to its first element!!!!"<<std::endl;
                    std::cout << p << std::endl;
                    p.resize(1);
                }
                else{
                    std::srand (std::time(NULL));// TODO: create either single seed for whole thing, or multiple seeds but keep track of them
                    //randomly choose a path
                    int r=std::rand()%possible_paths.size();
                    for (auto i=0;i<p.size();++i) p[i]=possible_paths[r][i];
                }
            }
        }
    }
    edge_migrations.close();

}