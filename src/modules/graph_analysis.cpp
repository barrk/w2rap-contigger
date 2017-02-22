//
// Created by Katie Barr (EI) on 22/02/2017.
//

#include "tclap/CmdLine.h"
#include "util/OutputLog.h"
#include "paths/HyperBasevector.h"


int main(const int argc, const char * argv[]) {

    TCLAP::CmdLine cmd("", ' ', "0.1");
    TCLAP::ValueArg<std::string> inFilename("i", "infile",
                                             "Name of input hbv", true, "", "string", cmd);
    TCLAP::ValueArg<std::string> outFilename("o", "outfile",
                                             "Name of output hbv", true, "", "string", cmd);

    std::string infile;
    std::string outfile;
    cmd.parse(argc, argv);
    // Get the value parsed by each arg.
    infile = inFilename.getValue();
    outfile = outFilename.getValue();

    HyperBasevector hbv;
    vec<int> hbvinv;

    OutputLog(2) <<"Loading graph..." << std::endl;
    BinaryReader::readFile(infile + ".hbv", &hbv);
    //Create inversion
    OutputLog(4) <<"Creating graph involution..." << std::endl;
    hbvinv.clear();
    hbv.Involution(hbvinv);
    // not sure what the difference between these 2 is - on ecoli small/large k, these give the same answer
    size_t n_components = hbv.NComponents();
    //int components = hbv.ConnectedComponents();

    vec< vec<int> > components;
    hbv.Components(components);
    OutputLog(1) << "Number of connected components: " << n_components << std::endl;

    std::cout << "Number of connected components: " << n_components << std::endl; //" components: " << components << std::endl;
    //auto self_loops = hbv.SelfLoops();

    //std::cout << "Number of self loops " << self_loops.size() << std::endl;

    vec<int> sources;
    vec<int> sinks;
    hbv.Sources(sources);
    hbv.Sinks(sinks);
    std::cout << "Number of sources: " << sources.size() << " and sinks " << sinks.size() << std::endl;

    // gather basic stats- number of nodes per component, average nodes per component, number of bases per component, average bases per component

    vec<size_t> sizes;
    vec<size_t> bases;
    int cycles = 0;
    size_t edge_bases_running_total = 0;
    size_t sizes_running_total = 0;
    size_t bases_running_total = 0;
    for (auto comp:components) {
        sizes.push_back(comp.size());
        for (auto edge: comp) {
            edge_bases_running_total += hbv.EdgeObject(edge).size();
        }
        if (hbv.HasCycle(comp)){
            cycles += 1;
        }
        bases.push_back(edge_bases_running_total);
        bases_running_total += edge_bases_running_total;
        sizes_running_total += comp.size();
        edge_bases_running_total = 0;
    }
    double mean_bases = bases_running_total/bases.size();
    double mean_sizes = sizes_running_total/sizes.size();

    // ideally want to classify components properly....
    std::cout << "number of components with cycles: " << cycles << std::endl;
    double bases_accum = 0.0;
    double sizes_accum = 0.0;
    std::pair< size_t, int > bases_max = std::make_pair(0,0);
    std::pair< size_t, int>  sizes_max = std::make_pair(0,0);
    for (int i=0; i < components.size(); i++){
        bases_accum += (bases[i] - mean_bases) * (bases[i] - mean_bases);
        sizes_accum += (sizes[i] - mean_sizes) * (sizes[i] - mean_sizes);
        bases_max = bases[i] > std::get<0>(bases_max) ? std::make_pair(bases[i], i) : bases_max;
        sizes_max = sizes[i] > std::get<0>(sizes_max) ? std::make_pair(bases[i], i) : sizes_max;
    }

    double bases_stdev = std::sqrt(bases_accum/(bases.size()-1));
    double sizes_stdev = std::sqrt(sizes_accum/(sizes.size()-1));

    std::cout << "Stats for bases per connected component, max: " << std::get<0>(bases_max) << " occurs at " << std::get<1>(bases_max) << ", mean: " << mean_bases << ", stdev: " <<  bases_stdev << std::endl;
    std::cout << "Stats for edges per connected component, max: " << std::get<0>(sizes_max) << " occurs at " << std::get<1>(sizes_max)<< ", mean: " << mean_sizes << ", stdev: " <<  sizes_stdev << std::endl;

    // the massive loopy part will be the part will be the component with the most edges, which should also be the component with the most bases, can now access this and focus analysis there..
    return 0;
}