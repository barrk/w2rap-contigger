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
    std::cout << "mapEdgesToLMPReads edge id to pair id map size lmp pairs size: " << edge_id_to_pair_id_map.size() << std::endl;
    return std::make_tuple(lmp_pairs, lmp_pairs_for_insert_size_estimation, edge_id_to_pair_id_map);
}


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
        auto edge_length = mHBV.EdgeObject(edge).size();
        // if we haven't traversed this edge, traverse it
        if (std::find(traversed_edge_list.begin(), traversed_edge_list.end(), edge) == traversed_edge_list.end()) {
            // if this edge takes us far enough away, add it to long fronteirs
            if (edge_length > (large_frontier_size - distance_traversed)){
                long_fronteirs.push_back(edge);
                // i want to store edges between start edge and edge, this stores attached edges, stores the ones traversed before twice for some readon
                traversed_edge_list.push_back(edge);
                std::vector<uint64_t> temp_path;
                for (auto intermediate_edge: intermediate_path){
                    temp_path.push_back(intermediate_edge);
                }
                paths_to_long_fronteirs.push_back(temp_path);
            } else { // if this edge does not take us far enough away, add its length to distance traversed go to next edge
                intermediate_path.push_back(edge);
                distance_traversed += mHBV.EdgeObject(edge).size();
                recursion_depth += 1;
                edges_beyond_distance(long_fronteirs, paths_to_long_fronteirs, intermediate_path, edge, traversed_edge_list, large_frontier_size,  recursion_depth, distance_traversed, direction);
                intermediate_path.clear();
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
    int same_in_out_degree = 0;
    int same_in_out_degree_complex = 0;
    int  solveable_regions_count = 0 ;
    std::map<uint64_t, int> mapping_counts;
    std::cout << "edge_id_to_pair_id_map size: " << edge_id_to_pair_id_map.size() << std::endl;
    std::vector<uint64_t > edges_to_look_for = {8, 321, 480, 391, 463, 358, 337, 478};
    std::cout << " " << std::endl;
    vector<uint64_t > traversed_edge_list;
    std::vector<std::vector<uint64_t> >  paths_to_long_fronteirs;
    std::vector<uint64_t> intermediate_path;
    std::vector<uint64_t>  long_frontiers_in;
    std::vector<uint64_t>  long_frontiers_out;
    std::vector<std::vector<uint64_t>> paths_to_separate;
    auto paths_containing_edges_to_look_for = {mEdgeToPathIds[8], mEdgeToPathIds[321], mEdgeToPathIds[480], mEdgeToPathIds[391],
                                                                 mEdgeToPathIds[463], mEdgeToPathIds[358], mEdgeToPathIds[337], mEdgeToPathIds[478]};

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

            if (long_frontiers_in.size() == long_frontiers_out.size()) {//} && long_frontiers_in.size() != 0) {
                same_in_out_degree += 1;
                for (auto path: paths_to_long_fronteirs){
                    std::cout << path_str(path) << std::endl;
                }
                if (long_frontiers_in.size() > 1) {
                    same_in_out_degree_complex += 1;
                    std::map<uint64_t, std::pair< std::vector<uint64_t>, std::vector<int> > > mapped_lmp_in;
                    std::map<uint64_t, std::pair< std::vector<uint64_t>, std::vector<int> > > mapped_lmp_out;
                    int count = 0;
                    // find number of lmp reads mapping to each large fronteir edge
                    for (auto edge: long_frontiers_in) {
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
                                    counts_for_each_solveabe_pair[std::make_pair(edge_in, edge_out)] += 1;

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

            }
            paths_to_long_fronteirs.clear();
            long_frontiers_in.clear();
            long_frontiers_out.clear();
    }
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