//
// Created by Katie Barr (EI) on 31/10/2016.
//

//
// Created by Bernardo Clavijo (TGAC) on 11/07/2016.
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

// actually just use dictionary from edge to path and pair ids
typedef struct {
    uint64_t edge;
    std::vector<uint64_t> intermediate_edges;
    std::vector<int> pair_ids;
} LMPmapping;

std::tuple <std::vector<LMPPair >, std::vector<LMPPair >, std::map<uint64_t, std::vector<int> > > PathFinderkb::mapEdgesToLMPReads(){
    KMatch kmatch(31);
    LMPMapper lmp_mapper(lmp_data, mHBV, mInv, kmatch);
    std::vector<LMPPair >  lmp_pairs;
    std::map<uint64_t, std::vector<int> > edge_id_to_pair_id_map;
    std::vector<LMPPair >  lmp_pairs_for_insert_size_estimation;
    lmp_mapper.LMPReads2MappedPairedEdgePaths(lmp_pairs, lmp_pairs_for_insert_size_estimation, edge_id_to_pair_id_map);
    std::cout << "mapEdgesToLMPReads lmp pairs size: " << lmp_pairs.size() << std::endl;
    std::cout << "mapEdgesToLMPReads edge id to pair id map size lmp pairs size: " << edge_id_to_pair_id_map.size() << std::endl;
    return std::make_tuple(lmp_pairs, lmp_pairs_for_insert_size_estimation, edge_id_to_pair_id_map);
}

// logic for standalone graph traversal
/*void traverse_graph(map <string, vector<tuple<string, string, string, bool> > > &edge_list, string start_node, int traversals,
            int total_traversals, ofstream &outfile, vector<pair<string, string>> &traversed_edge_list){
            traversals += 1;
            vector<tuple<string, string, string, bool> > adjacent_nodes = edge_list[start_node];
            cout << adjacent_nodes.size();
            for (auto node = adjacent_nodes.begin(); node != adjacent_nodes.end(); ++node){ //= adjacent_nodes.begin(); node != adjacent_nodes.end(); ++node){
                string end_node = get<0>(*node);
                pair<string, string> edge_set = make_pair(start_node, end_node);
                if (find(traversed_edge_list.begin(), traversed_edge_list.end(), edge_set) == traversed_edge_list.end()){
                    //determine whether vertex order must be flipped to ensure sequences concatenated in correct order
                    bool edge_direction_correct = get<3>(*node);
                    if (edge_direction_correct){
                        outfile << "L" << "\t" << start_node << "\t" << get<1>(*node)  << "\t" << end_node << "\t" << get<2>(*node) << "\t" << "0M" << endl;
                        } else {
                        outfile << "L" << "\t" << end_node << "\t" << get<2>(*node)  << "\t" << start_node << "\t" << get<1>(*node) << "\t" << "0M" << endl;
                        }
                    traversed_edge_list.push_back(edge_set);
                    if ((traversals+1) <= total_traversals){
                        traverse_graph(edge_list, get<0>(*node), traversals, total_traversals, outfile, traversed_edge_list);
                }
            }
}*/


void PathFinderkb::edges_beyond_distance(std::vector<uint64_t>  & long_fronteirs, std::vector<std::vector<uint64_t> >  & paths_to_long_fronteirs, std::vector<uint64_t> & intermediate_path, uint64_t e, std::vector<uint64_t > & traversed_edge_list, uint64_t large_frontier_size, int recursion_depth=0, int distance_traversed=0, std::string direction="right") {
    /*
     * to find edges in region lmp pairs might solve, we want to traverse graph until we're large_frontier_size away from e
     */
    std::vector<uint64_t> edges;
    if (direction=="right"){edges = prev_edges[e];
    } else {
        edges = next_edges[e];
    }
    //std::cout << "traversal called for edge: " << e << std::endl;
    // edges between the start edge and the end edge are ones we traverse through in that recursion stack,
    // only want paths from beginning node, i.e. when recursion depth is 0, but paths can split at any time
    for (auto edge:edges){
        /*if (recursion_depth==0){
            std::vector<uint64_t> intermediate_path;
            paths_to_long_fronteirs.push_back(intermediate_path); this doesn't enbla eus to ad paths which fork later
        }*/
        //std::cout << "traversing edge: " << edge << " distance traverse: " << distance_traversed <<std::endl;
        auto edge_length = mHBV.EdgeObject(edge).size();
        // if we haven't traversed this edge, traverse it
        if (std::find(traversed_edge_list.begin(), traversed_edge_list.end(), edge) == traversed_edge_list.end()) {
            //std::cout << "edge not traversed, length: " << edge_length << " dustance traversed "  << distance_traversed << " distance left to traverse: " << (large_frontier_size - distance_traversed) <<std::endl;
            // if this edge takes us far enough away, add it to long fronteirs
            if (edge_length > (large_frontier_size - distance_traversed)){
                long_fronteirs.push_back(edge);
                // i want to store edges between start edge and edge, this stores attached edges, stores the ones traversed before twice for some readon
                traversed_edge_list.push_back(edge);
                //std::cout << "fonteir found, edge: " << edge << std::endl;
                std::vector<uint64_t> temp_path;
                for (auto intermediate_edge: intermediate_path){
                    temp_path.push_back(intermediate_edge);
                }
                paths_to_long_fronteirs.push_back(temp_path);
                //std::cout << "intermediate edges:" << path_str(paths_to_long_fronteirs.back()) << std::endl;
                //std::cout << path_str(paths_to_long_fronteirs.back()) << std::endl;
            } else { // if this edge does not take us far enough away, add its length to distance traversed go to next edge
                //std::cout << "going to next edge"  << std::endl;
                intermediate_path.push_back(edge);
                distance_traversed += mHBV.EdgeObject(edge).size();
                recursion_depth += 1;
                edges_beyond_distance(long_fronteirs, paths_to_long_fronteirs, intermediate_path, edge, traversed_edge_list, large_frontier_size,  recursion_depth, distance_traversed, direction);
                intermediate_path.clear();
            }
        }
    }
    //return long_fronteirs;

}

void PathFinderkb::traverseGraphCreateSubgraph(vector<uint64_t>  edges_to_add, std::map<int, uint64_t> & vertex_subgraph_map, HyperBasevector & subgraph,
                                      vector<uint64_t > & traversed_edge_list, int edges_to_traverse, int traversals=0) {
    for (auto edge = edges_to_add.begin();
         edge != edges_to_add.end(); ++edge) { //= adjacent_nodes.begin(); node != adjacent_nodes.end(); ++node){
        if (std::find(traversed_edge_list.begin(), traversed_edge_list.end(), *edge) == traversed_edge_list.end()) {
            basevector bv = mHBV.EdgeObject(*edge);
            int prev_vertex = subgraph.N();
            int next_vertex = subgraph.N() + 1;
            vertex_subgraph_map[prev_vertex] = mToLeft[*edge];
            vertex_subgraph_map[next_vertex] = mToRight[*edge];
            subgraph.AddVertices(2);
            subgraph.AddEdge(prev_vertex, next_vertex, bv);
            traversed_edge_list.push_back(*edge);
            if ((traversals + 1) <= edges_to_traverse) {
                traverseGraphCreateSubgraph(edges_to_add, vertex_subgraph_map, subgraph, traversed_edge_list,
                                            edges_to_traverse,
                                            traversals);
            }
        }
    }

}


void PathFinderkb::addEdgeToSubGraph(int edge_to_add, std::map<int, uint64_t> & vertex_subgraph_map, HyperBasevector & subgraph, vector<uint64_t > & traversed_edge_list){
    // for each edge to add, need to add vertices and n adjacent edges
    auto bv = mHBV.EdgeObject(edge_to_add);
    // store relationship between vertices here and in main graph
    int prev_vertex = subgraph.N();
    int next_vertex = subgraph.N() + 1;
    vertex_subgraph_map[prev_vertex] = mToLeft[edge_to_add];
    vertex_subgraph_map[next_vertex] = mToRight[edge_to_add];
    subgraph.AddVertices(2);
    subgraph.AddEdge(prev_vertex,  next_vertex, bv);// add edge takes vertex from, vertex to and base vector
    vector<uint64_t> edges_to_add;
    edges_to_add.reserve(prev_edges[edge_to_add].size() + next_edges[edge_to_add].size());
    edges_to_add.insert(edges_to_add.end(), prev_edges[edge_to_add].begin(), prev_edges[edge_to_add].end());
    edges_to_add.insert(edges_to_add.end(), next_edges[edge_to_add].begin(), next_edges[edge_to_add].end());
    // edges_to_add, std::map<int, uint64_t> & vertex_subgraph_map, HyperBasevector & subgraph,
    //vector<uint64_t > & traversed_edge_list, int edges_to_traverse, int traversals=0
    traverseGraphCreateSubgraph(edges_to_add, vertex_subgraph_map,  subgraph,
                                traversed_edge_list,  4,  0);
}


void PathFinderkb::gatherStats() {

    int large_frontier_size = 3500;
    init_prev_next_vectors();// mToLeft and mToright are to/from vertices
    std::tuple<std::vector<LMPPair >, std::vector<LMPPair >, std::map<uint64_t, std::vector<int> > > mapping_results = mapEdgesToLMPReads();
    std::vector<LMPPair > pairs_for_scaffolding = std::get<0>(mapping_results);
    std::vector<LMPPair > pairs_for_insert_size_calculation = std::get<1>(mapping_results);
    std::map<uint64_t, std::vector<int> >  edge_id_to_pair_id_map = std::get<2>(mapping_results);
    // can use linux tools to get number of reads mappig to these edges to save time searching in here
    // once we have the read ids, can select a suitable subset to run quickly during development
    std::ofstream edge_ids_with_long_frontiera;
    edge_ids_with_long_frontiera.open("/Users/barrk/Documents/arabidopsis_data/ecoli_fronteir_edge_is_with_same_in_out_degree_single_edge_mappings.txt");
    std::ofstream solveable_regions;
    solveable_regions.open("/Users/barrk/Documents/ecoli_data/solveable_regions_single_edge_mappings.txt");
    std::ofstream paths_over_subgraph;
    paths_over_subgraph.open("/Users/barrk/Documents/ecoli_data/solveable_regions_subgraph_single_edge_mappings.txt");
    int same_in_out_degree = 0;
    int same_in_out_degree_complex = 0;
    int  solveable_regions_count = 0 ;
    std::map<uint64_t, int> mapping_counts;
    std::cout << "edge_id_to_pair_id_map size: " << edge_id_to_pair_id_map.size() << std::endl;
    std::vector<uint64_t > edges_to_look_for = {8, 321, 480, 391, 463, 358, 337, 478};
    /*for (auto edge: edges_to_look_for){
        std::cout << "pair ids for edge: " << edge << std::endl;
        for (auto pid: edge_id_to_pair_id_map[edge]){
            std::cout << pid << " " << std::endl;
        }

    }*/
    std::cout << " " << std::endl;
    vector<uint64_t > traversed_edge_list;
    std::vector<std::vector<uint64_t> >  paths_to_long_fronteirs;
    std::vector<uint64_t> intermediate_path;
    std::vector<uint64_t>  long_frontiers_in;
    std::vector<uint64_t>  long_frontiers_out;
    std::vector<std::vector<uint64_t>> paths_to_separate;
    auto paths_containing_edges_to_look_for = {mEdgeToPathIds[8], mEdgeToPathIds[321], mEdgeToPathIds[480], mEdgeToPathIds[391],
                                                                 mEdgeToPathIds[463], mEdgeToPathIds[358], mEdgeToPathIds[337], mEdgeToPathIds[478]};
    /*for (auto edge:edges_to_look_for) {
        std::cout << "edge: " << edge << std::endl;
        for (auto path_id: mEdgeToPathIds[edge]) {
            // when printing path id here, it printed absolutely loads
            if (mPaths[path_id].size() > 2) {
                std::cout << "path id: " << path_id << std::endl;
                std::vector<uint64_t> pathvec;
                for (auto e:mPaths[path_id]) { pathvec.push_back(e); }
                std::cout << path_str(pathvec) << std::endl;
            }
        }
    }*/
    for (auto edge_index: edges_to_look_for) {
        if (std::find(edges_to_look_for.begin(), edges_to_look_for.end(), edge_index) != edges_to_look_for.end()) {
            std::cout << "edge of interest:" << edge_index << std::endl;
        }
        // e < mInv[e] checks that this is a forward directed edge? nope, canonical representation
        //if (edge_index < mInv[edge_index]) {//&& mHBV.EdgeObject(edge_index).size() < large_frontier_size) {
            edges_beyond_distance(long_frontiers_in, paths_to_long_fronteirs, intermediate_path, edge_index, traversed_edge_list, large_frontier_size, 0, 0, "right");
            traversed_edge_list.clear();
            std::cout << "traversed edge list cleared: " << traversed_edge_list.size() << std::endl;
            edges_beyond_distance(long_frontiers_out, paths_to_long_fronteirs, intermediate_path, edge_index, traversed_edge_list, large_frontier_size, 0, 0, "left");

            if (edge_index == 321) {
                std::cout << "long frontiers of edge 321" << std::endl;
                for (auto e: long_frontiers_in) {
                    std::cout << e << " ";
                }
                std::cout << std::endl;
                for (auto e: long_frontiers_out) {
                    std::cout << e << " ";
                }
                std::cout << std::endl;

            }
            if (long_frontiers_in.size() == long_frontiers_out.size()) {//} && long_frontiers_in.size() != 0) {
                same_in_out_degree += 1;
                edge_ids_with_long_frontiera << "Edge: " << edge_index << " size: " << long_frontiers_in.size()
                                             << std::endl;
                for (auto path: paths_to_long_fronteirs){
                    std::cout << path_str(path) << std::endl;
                }
                if (long_frontiers_in.size() > 1) {
                    std::cout << "looking for mappings to edges in fronteir: " << std::endl;
                    for (auto frontier: long_frontiers_in) {
                        std::cout << frontier << " ";
                    }
                    std::cout << std::endl;
                    std::cout << "to out fronteir: " << std::endl;
                    for (auto frontier: long_frontiers_out) {
                        std::cout << frontier << " ";
                    }
                    std::cout << std::endl;
                    same_in_out_degree_complex += 1;
                    std::map<uint64_t, std::pair< std::vector<uint64_t>, std::vector<int> > > mapped_lmp_in;
                    std::map<uint64_t, std::pair< std::vector<uint64_t>, std::vector<int> > > mapped_lmp_out;
                    int count = 0;
                    // find number of lmp reads mapping to each large fronteir edge
                    for (auto edge: long_frontiers_in) {
                        std::cout << "Edge in: " <<edge <<std::endl;
                        std::cout <<  path_str(paths_to_long_fronteirs[count]) <<std::endl;
                        std::vector<uint64_t> intermediate_edges;
                        for (auto int_edge: paths_to_long_fronteirs[count]){
                            intermediate_edges.push_back(int_edge);
                        }
                        std::vector<int> pair_ids;
                        for (auto pair_id :  edge_id_to_pair_id_map[edge]) {
                            //std::cout  << pair_id << " ";
                            pair_ids.push_back(pair_id);
                            mapping_counts[edge] += 1;
                        }
                        count += 1;
                        std::pair< std::vector<uint64_t>, std::vector<int> > result = std::make_pair(intermediate_edges, pair_ids);
                        mapped_lmp_in[edge] = result;
                    }
                    // if there is
                    for (auto edge: long_frontiers_out) {
                        std::cout << "Edge out: " <<edge <<std::endl;
                        std::cout <<  path_str(paths_to_long_fronteirs[count]) <<std::endl;
                        std::vector<uint64_t> intermediate_edges;
                        for (auto int_edge: paths_to_long_fronteirs[count]){
                            intermediate_edges.push_back(int_edge);
                        }
                        std::vector<int> pair_ids;
                        for (auto pair_id :  edge_id_to_pair_id_map[mInv[edge]]) {
                            //std::cout  << pair_id << " ";
                            mapping_counts[edge] += 1;
                            pair_ids.push_back(pair_id);

                        }
                        std::cout << std::endl;
                        count += 1;
                        std::pair< std::vector<uint64_t>, std::vector<int> > result = std::make_pair(intermediate_edges, pair_ids);
                        mapped_lmp_out[edge] = result;

                    }
                    /*if (edge_index == 321) {
                        std::cout << "mapping out/in, if there is a pair id in common, that pair resolves a path"
                                  << std::endl;
                        for (auto lmp_pair_in: mapped_lmp_in) {
                            std::cout << lmp_pair_in << " ";
                        }
                        std::cout << "pair ids mapping to edge leavig complex regions" << std::endl;
                        for (auto lmp_pair_out: mapped_lmp_out) {
                            std::cout << lmp_pair_out << " ";
                        }

                    }*/
                    std::map<std::pair<uint64_t, uint64_t>, int > counts_for_each_solveabe_pair;
                    // find pairs which are in both mapped_lmp_in, and mapped_lmp_out- these can be used to solve this region
                    for (auto edge_in: long_frontiers_in) {
                        // TODO this way is horrible, need to design a better data structure
                        vector<int> p1_ids = mapped_lmp_in[edge_in].second;
                        for (auto p1_id: p1_ids){
                            for (auto edge_out: long_frontiers_out) {
                                vector<int> p2_ids = mapped_lmp_out[edge_out].second;
                                if ((std::find(p2_ids.begin(), p2_ids.end(), p1_id) !=
                                        p2_ids.end())) {
                                    //std::cout << "pair id: " << p1_id << std::endl;
                                    //auto pair = pairs_for_scaffolding[p1_id];
                                    counts_for_each_solveabe_pair[std::make_pair(edge_in, edge_out)] += 1;
                                    /*for (auto p:pair.p1) {
                                        std::cout << p << " ";
                                    }
                                    for (auto p:pair.p2) {
                                        std::cout << p << " ";
                                    }
                                    std::cout << std::endl;*/

                                }
                            }
                        }
                    }
                    for (auto edge_pair_count: counts_for_each_solveabe_pair){
                        std::pair<uint64_t, uint64_t> key = edge_pair_count.first;
                        uint64_t edge_in = key.first;
                        uint64_t edge_out = key.second;
                        int count = edge_pair_count.second;
                        std::cout << "count for pair:" << edge_pair_count.first.first <<  " " << edge_pair_count.first.second << " : " <<count <<std::endl;
                        if (count > 5){
                            std::vector<uint64_t> path_in = mapped_lmp_in[edge_in].first;
                            std::vector<uint64_t> path_out = mapped_lmp_out[edge_out].first;
                            std::vector<uint64_t> full_path;
                            full_path.push_back(edge_in);
                            std::cout << "path in " << path_str(path_in) << std::endl;
                            std::cout << "path out " << path_str(path_out) << std::endl;
                            for (auto e:path_in){
                                full_path.push_back(e);
                            }
                            full_path.push_back(edge_index);
                            for (auto e:path_out){
                                full_path.push_back(e);
                            }
                            full_path.push_back(edge_out);
                            std::cout << "path to separate: " << path_str(full_path) << std::endl;
                            paths_to_separate.push_back(full_path);
                        }
                    }
                }
                //}

            }
            paths_to_long_fronteirs.clear();
            long_frontiers_in.clear();
            long_frontiers_out.clear();
        //}
    }
    edge_ids_with_long_frontiera.close();
    solveable_regions.close();
    paths_over_subgraph.close();
    std::cout << "Regions with same degree in and out:" << same_in_out_degree << std::endl;
    std::cout << "Complex regions with same degree in and out:" << same_in_out_degree_complex << std::endl;
    std::cout << "Complex regions solveable wit lmp reads:" << solveable_regions_count<< std::endl;
    uint64_t sep=0;
    std::map<uint64_t,std::vector<uint64_t>> old_edges_to_new;
    for (auto p:paths_to_separate){

        if (old_edges_to_new.count(p.front()) > 0 or old_edges_to_new.count(p.back()) > 0) {
            std::cout<<"WARNING: path starts or ends in an already modified edge, skipping"<<std::endl;
            continue;
        }

        // this creates new vertices at the end of the graph, connects start to these, and removes old connection then repeats on the involution
        auto oen=separate_path(p, true);
        if (oen.size()>0) {
            for (auto et:oen){
                if (old_edges_to_new.count(et.first)==0) old_edges_to_new[et.first]={};
                for (auto ne:et.second) old_edges_to_new[et.first].push_back(ne);
            }
            sep++;
        }
    }
    if (old_edges_to_new.size()>0) {
        migrate_readpaths(old_edges_to_new);
    }
    std::cout<<" "<<sep<<" paths separated!"<<std::endl;
    BinaryWriter::writeFile("/Users/barrk/Documents/ecoli_dataset/v1/all_paths_separated.hbv", mHBV);
    WriteReadPathVec(mPaths,"/Users/barrk/Documents/ecoli_dataset/v1/all__paths_separated.paths");
    std::cout << "Dumping gfa" << std::endl;
    int MAX_CELL_PATHS = 50;
    int MAX_DEPTH = 10;
    GFADump("/Users/barrk/Documents/ecoli_dataset/v1/all_paths_separated", mHBV, mInv, mPaths, MAX_CELL_PATHS, MAX_DEPTH, true);


}

void PathFinderkb::resolveComplexRegionsUsingLMPData(){
    std::tuple<std::vector<LMPPair >, std::vector<LMPPair >, std::map<uint64_t, std::vector<int> > > mapping_results = mapEdgesToLMPReads();
    std::vector<LMPPair > lmp_pairs = std::get<0>(mapping_results);
    std::vector<LMPPair > pairs_for_insert_size_calculation = std::get<1>(mapping_results);
    std::map<uint64_t, std::vector<int> >  edge_id_to_pair_id_map = std::get<2>(mapping_results);
    std::cout << "lmp pairs size: " << lmp_pairs.size() << std::endl;
    // then find complex regions, for now in same way as original pathfinder
    init_prev_next_vectors();
    std::cout<<"vectors initialised"<<std::endl;
    std::vector<std::vector<uint64_t>> paths_to_separate;
    std::set<int> paths_found;
    bool reversed = false;
    // simplest approach is to separate all paths from LMP pairs
    for (auto pair: lmp_pairs) {
        auto p1 = pair.p1;
        auto p2 = pair.p2;
        for (auto e1: p1) {// these are the edges
            for (auto e2: p2) {
                // p1 and p2 might be mapped with either first
                // the second part of this conditional vaguely determines whether either of the edges is complex
                if ((e1 != e2) && (next_edges[e1].size() > 1 || prev_edges[e2].size() > 1 || next_edges[e2].size() > 1 ||
                    prev_edges[e1].size() > 1)) {
                    //std::cout << "pair spans complex region" << std::endl;
                    auto shared_paths = 0;
                    for (auto inp:mEdgeToPathIds[e1]) {// find paths associated with in_e
                        //std::cout << "inp" << inp << std::endl;
                        for (auto outp:mEdgeToPathIds[e2]) {// ditto out edge
                            //std::cout << "outp" << outp << std::endl;
                            if ((inp == outp) && (paths_found.find(inp) ==
                                                  paths_found.end())) {// if they're on the same path and this path hasn't been added
                                paths_found.insert(inp);
                                std::cout << "shared path found" << std::endl;
                                shared_paths++;
                                if (shared_paths == 1) {
                                    // need to determine which of e1 and e2 is first i npath for indexes below, mPaths doesn't have find method on it, even though really its a vector
                                    bool seen_e1 = false;
                                    int16_t i = 0;
                                    while (mPaths[inp][i] != e2) {
                                        if (mPaths[inp][i] == e1){ seen_e1 = true;}
                                        i++;
                                    }
                                    auto first = seen_e1 ? e1:e2;
                                    auto second = seen_e1 ? e2:e1;
                                    std::vector<uint64_t> pv;
                                    for (auto e:mPaths[inp]) pv.push_back(e); // create vector of edges on theshared pat
                                    std::cout << "found first path from " << e1 << " to " << e2 << path_str(pv)
                                              << std::endl;
                                    paths_to_separate.push_back({});
                                    int16_t ei = 0;
                                    while (mPaths[inp][ei] != first) ei++; // find which edge on the path is the in edge
                                    // add all edges on the path until the out edge
                                    while (mPaths[inp][ei] != second && ei < mPaths[inp].size())
                                        paths_to_separate.back().push_back(mPaths[inp][ei++]);
                                    if (ei >= mPaths[inp].size()) {
                                        std::cout << "reversed path detected!" << std::endl;
                                        reversed = true;
                                    }
                                    paths_to_separate.back().push_back(second);
                                    std::cout<<"added!"<< path_str(paths_to_separate.back()) << std::endl;
                                }
                            }
                            if (shared_paths > 1) {
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

        uint64_t sep=0;
        std::map<uint64_t,std::vector<uint64_t>> old_edges_to_new;
        for (auto p:paths_to_separate){

            if (old_edges_to_new.count(p.front()) > 0 or old_edges_to_new.count(p.back()) > 0) {
                std::cout<<"WARNING: path starts or ends in an already modified edge, skipping"<<std::endl;
                continue;
            }
            std::cout<<"separating path  "<<path_str(p)<< std::endl;

            auto oen=separate_path(p, true);
            if (oen.size()>0) {
                for (auto et:oen){
                    if (old_edges_to_new.count(et.first)==0) old_edges_to_new[et.first]={};
                    for (auto ne:et.second) old_edges_to_new[et.first].push_back(ne);
                }
                sep++;
            }
        }
        if (old_edges_to_new.size()>0) {
            migrate_readpaths(old_edges_to_new);
        }
        std::cout<<" "<<sep<<" paths separated!"<<std::endl;
    BinaryWriter::writeFile("/Users/barrk/Documents/ecoli_dataset/after_lmp_mapping.hbv", mHBV);
    WriteReadPathVec(mPaths, "/Users/barrk/Documents/ecoli_dataset/after_lmp_mapping.paths");

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


void PathFinderkb::untangle_complex_in_out_choices(uint64_t large_frontier_size, bool verbose_separation) {
    //find a complex path
    uint64_t qsf=0,qsf_paths=0;
    uint64_t msf=0,msf_paths=0;
    init_prev_next_vectors();// need to clarify what mToLeft and mToRight are, guessing they're just into/out of each node
    std::cout<<"vectors initialised"<<std::endl;
    std::set<std::array<std::vector<uint64_t>,2>> seen_frontiers,solved_frontiers;
    std::vector<std::vector<uint64_t>> paths_to_separate;
    for (int e = 0; e < mHBV.EdgeObjectCount(); ++e) {
        if (e < mInv[e] && mHBV.EdgeObject(e).size() < large_frontier_size) {
            auto f=get_all_long_frontiers(e, large_frontier_size);// return edges in and out of this edge which represent long sequences of bases
            // think f[0] is vec of edge ids going into this edge, and f[1] is vec of edge ids going out
            if (f[0].size()>1 and f[1].size()>1 and seen_frontiers.count(f)==0){ // this determines whether it s acomplex edge
                // if theres more than 1 edge in, and more than 1 edge out, and we haven't already dealt with this in/out combination
                seen_frontiers.insert(f);
                bool single_dir=true;
                // this determines whether the same edge is both going into and out of current edge- i.e. its a palindrome
                for (auto in_e:f[0]) for (auto out_e:f[1]) if (in_e==out_e) {single_dir=false;break;}
                if (single_dir) {
                    std::cout<<" Single direction frontiers for complex region on edge "<<e<<" IN:"<<path_str(f[0])<<" OUT: "<<path_str(f[1])<<std::endl;
                    std::vector<int> in_used(f[0].size(),0);// initialise vec for each in edge
                    std::vector<int> out_used(f[1].size(),0); //  ditto out edge
                    std::vector<std::vector<uint64_t>> first_full_paths;
                    bool reversed=false;
                    for (auto in_i=0;in_i<f[0].size();++in_i) {
                        auto in_e=f[0][in_i];
                        for (auto out_i=0;out_i<f[1].size();++out_i) {
                            auto out_e=f[1][out_i];
                            auto shared_paths = 0;
                            for (auto inp:mEdgeToPathIds[in_e]) // find paths associated with in_e
                                for (auto outp:mEdgeToPathIds[out_e]) // ditto out edge
                                    if (inp == outp) {// if they're on the same path

                                        shared_paths++;
                                        if (shared_paths==1){
                                            std::vector<uint64_t > pv;
                                            for (auto e:mPaths[inp]) pv.push_back(e); // create vector of edges on theshared pat
                                            std::cout<<"found first path from "<<in_e<<" to "<< out_e << path_str(pv)<< std::endl;
                                            first_full_paths.push_back({});
                                            int16_t ei=0;
                                            while (mPaths[inp][ei]!=in_e) ei++; // find which edge on the path is the in edge
                                            // add all edges on the path until the out edge
                                            while (mPaths[inp][ei]!=out_e && ei<mPaths[inp].size()) first_full_paths.back().push_back(mPaths[inp][ei++]);
                                            if (ei>=mPaths[inp].size()) {
                                                std::cout<<"reversed path detected!"<<std::endl;
                                                reversed=true;
                                            }
                                            first_full_paths.back().push_back(out_e);
                                            //std::cout<<"added!"<<std::endl;
                                        }
                                    }
                            //check for reverse paths too- same but on inovlution
                            for (auto inp:mEdgeToPathIds[mInv[out_e]])
                                for (auto outp:mEdgeToPathIds[mInv[in_e]])
                                    if (inp == outp) {

                                        shared_paths++;
                                        if (shared_paths==1){//not the best solution, but should work-ish
                                            std::vector<uint64_t > pv;
                                            for (auto e=mPaths[inp].rbegin();e!=mPaths[inp].rend();++e) pv.push_back(mInv[*e]);
                                            std::cout<<"found first path from "<<in_e<<" to "<< out_e << path_str(pv)<< std::endl;
                                            first_full_paths.push_back({});
                                            int16_t ei=0;
                                            while (pv[ei]!=in_e) ei++;

                                            while (pv[ei]!=out_e && ei<pv.size()) first_full_paths.back().push_back(pv[ei++]);
                                            if (ei>=pv.size()) {
                                                std::cout<<"reversed path detected!"<<std::endl;
                                                reversed=true;
                                            }
                                            first_full_paths.back().push_back(out_e);
                                            //std::cout<<"added!"<<std::endl;
                                        }
                                    }
                            if (shared_paths) {
                                out_used[out_i]++;
                                in_used[in_i]++;
                                //std::cout << "  Shared paths " << in_e << " --> " << out_e << ": " << shared_paths << std::endl;

                            }
                        }
                    }
                    if ((not reversed) and std::count(in_used.begin(),in_used.end(),1) == in_used.size() and
                        std::count(out_used.begin(),out_used.end(),1) == out_used.size()){
                        std::cout<<" REGION COMPLETELY SOLVED BY PATHS!!!"<<std::endl;
                        solved_frontiers.insert(f);
                        for (auto p:first_full_paths) paths_to_separate.push_back(p);
                    }
                }

            }
        }
    }
    std::cout<<"Complex Regions solved by paths: "<<solved_frontiers.size() <<"/"<<seen_frontiers.size()<<" comprising "<<paths_to_separate.size()<<" paths to separate"<< std::endl;
    //std::cout<<"Complex Regions quasi-solved by paths (not acted on): "<< qsf <<"/"<<seen_frontiers.size()<<" comprising "<<qsf_paths<<" paths to separate"<< std::endl;
    //std::cout<<"Multiple Solution Regions (not acted on): "<< msf <<"/"<<seen_frontiers.size()<<" comprising "<<msf_paths<<" paths to separate"<< std::endl;

    uint64_t sep=0;
    std::map<uint64_t,std::vector<uint64_t>> old_edges_to_new;
    for (auto p:paths_to_separate){

        if (old_edges_to_new.count(p.front()) > 0 or old_edges_to_new.count(p.back()) > 0) {
            std::cout<<"WARNING: path starts or ends in an already modified edge, skipping"<<std::endl;
            continue;
        }

        // this creates new vertices at the end of the graph, connects start to these, and removes old connection then repeats on the involution
        auto oen=separate_path(p, verbose_separation);
        if (oen.size()>0) {
            for (auto et:oen){
                if (old_edges_to_new.count(et.first)==0) old_edges_to_new[et.first]={};
                for (auto ne:et.second) old_edges_to_new[et.first].push_back(ne);
            }
            sep++;
        }
    }
    if (old_edges_to_new.size()>0) {
        migrate_readpaths(old_edges_to_new);
    }
    std::cout<<" "<<sep<<" paths separated!"<<std::endl;
}



std::array<std::vector<uint64_t>,2> PathFinderkb::get_all_long_frontiers(uint64_t e, uint64_t large_frontier_size){
    //TODO: return all components in the region
    std::set<uint64_t> seen_edges, to_explore={e}, in_frontiers, out_frontiers;

    while (to_explore.size()>0){
        std::set<uint64_t> next_to_explore;
        for (auto x:to_explore){ //to_explore: replace rather and "update" (use next_to_explore)

            if (!seen_edges.count(x)){

                //What about reverse complements and paths that include loops that "reverse the flow"?
                if (seen_edges.count(mInv[x])) return std::array<std::vector<uint64_t>,2>(); //just cancel for now

                for (auto p:prev_edges[x]) {
                    if (mHBV.EdgeObject(p).size() >= large_frontier_size )  {
                        //What about frontiers on both sides?
                        in_frontiers.insert(p);
                        for (auto other_n:next_edges[p]){
                            if (!seen_edges.count(other_n)) {
                                if (mHBV.EdgeObject(other_n).size() >= large_frontier_size) {
                                    out_frontiers.insert(other_n);
                                    seen_edges.insert(other_n);
                                }
                                else next_to_explore.insert(other_n);
                            }
                        }
                    }
                    else if (!seen_edges.count(p)) next_to_explore.insert(p);
                }

                for (auto n:next_edges[x]) {
                    if (mHBV.EdgeObject(n).size() >= large_frontier_size) {
                        //What about frontiers on both sides?
                        out_frontiers.insert(n);
                        for (auto other_p:prev_edges[n]){
                            if (!seen_edges.count(other_p)) {
                                if (mHBV.EdgeObject(other_p).size() >= large_frontier_size) {
                                    in_frontiers.insert(other_p);
                                    seen_edges.insert(other_p);
                                }
                                else next_to_explore.insert(other_p);
                            }
                        }
                    }
                    else if (!seen_edges.count(n)) next_to_explore.insert(n);
                }
                seen_edges.insert(x);
            }
            if (seen_edges.size()>50) {
                return std::array<std::vector<uint64_t>,2>();
            }

        }
        to_explore=next_to_explore;
    }

    if (to_explore.size()>0) return std::array<std::vector<uint64_t>,2>();
    //the "canonical" representation is the one that has the smalled edge on the first vector, and bot ordered

    if (in_frontiers.size()>0 and out_frontiers.size()>0) {
        uint64_t min_in=*in_frontiers.begin();
        for (auto i:in_frontiers){
            if (i<min_in) min_in=i;
            if (mInv[i]<min_in) min_in=mInv[i];
        }
        uint64_t min_out=*out_frontiers.begin();
        for (auto i:out_frontiers){
            if (i<min_out) min_out=i;
            if (mInv[i]<min_out) min_out=mInv[i];
        }
        if (min_out<min_in){
            std::set<uint64_t> new_in, new_out;
            for (auto e:in_frontiers) new_out.insert(mInv[e]);
            for (auto e:out_frontiers) new_in.insert(mInv[e]);
            in_frontiers=new_in;
            out_frontiers=new_out;
        }
    }

    //std::sort(in_frontiers.begin(),in_frontiers.end());
    //std::sort(out_frontiers.begin(),out_frontiers.end());
    std::array<std::vector<uint64_t>,2> frontiers={std::vector<uint64_t>(in_frontiers.begin(),in_frontiers.end()),std::vector<uint64_t>(out_frontiers.begin(),out_frontiers.end())};



    return frontiers;
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
    if (verbose_separation) std::cout<<"Migrating edge "<<mInv[p[0]]<<" From node old: "<<mToLeft[mInv[p[0]]]<<" new: "<<current_vertex_rev<<std::endl;
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
        if (verbose_separation)  std::cout<<"Edge "<<nef<<": copy of "<<p[ei]<<": "<<prev_vertex_fw<<" - "<<current_vertex_fw<<std::endl;
        mToLeft.push_back(prev_vertex_fw);
        mToRight.push_back(current_vertex_fw);
        if (! old_edges_to_new.count(p[ei]))  old_edges_to_new[p[ei]]={};
        old_edges_to_new[p[ei]].push_back(nef);

        auto ner=mHBV.AddEdge(current_vertex_rev,prev_vertex_rev,mHBV.EdgeObject(mInv[p[ei]]));
        if (verbose_separation) std::cout<<"Edge "<<ner<<": copy of "<<mInv[p[ei]]<<": "<<current_vertex_rev<<" - "<<prev_vertex_rev<<std::endl;
        mToLeft.push_back(current_vertex_rev);
        mToRight.push_back(prev_vertex_rev);
        if (! old_edges_to_new.count(mInv[p[ei]]))  old_edges_to_new[mInv[p[ei]]]={};
        old_edges_to_new[mInv[p[ei]]].push_back(ner);

        mInv.push_back(ner);
        mInv.push_back(nef);
        mEdgeToPathIds.resize(mEdgeToPathIds.size()+2);
    }
    if (verbose_separation) std::cout<<"Migrating edge "<<p[p.size()-1]<<" From node old: "<<mToLeft[p[p.size()-1]]<<" new: "<<current_vertex_fw<<std::endl;
    mHBV.GiveEdgeNewFromVx(p[p.size()-1],mToLeft[p[p.size()-1]],current_vertex_fw);// attach new edges back into graph at right position
    //mToLeft[p[p.size()-1]] = current_vertex_fw;
    if (verbose_separation) std::cout<<"Migrating edge "<<mInv[p[p.size()-1]]<<" To node old: "<<mToRight[mInv[p[p.size()-1]]]<<" new: "<<current_vertex_rev<<std::endl;
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
                    std::cout << "Migrating edge: " << p[i] << " to edge " << possible_new_edges[i][0] << std::endl;
                    edge_migrations << "Migrating edge: " << p[i] << " to edge " << possible_new_edges[i][0] << std::endl;
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
                    p.resize(1);
                }
                else{
                    std::srand (std::time(NULL));
                    //randomly choose a path
                    int r=std::rand()%possible_paths.size();
                    for (auto i=0;i<p.size();++i) p[i]=possible_paths[r][i];
                }
            }
        }
    }
    edge_migrations.close();

}