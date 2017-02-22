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
    size_t n_components = hbv.NComponents();
    OutputLog(1) << "Number of connected components: " << n_components << std::endl;

    std::cout << "Number of connected components: " << n_components << std::endl;

    
    return 0;
}