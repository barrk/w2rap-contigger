//
// Created by Katie Barr (EI) and Rob Vickerstaff (EMR) on 24/02/2017.
//

#include "paths/long/large/ExtractReads.h"
#include "tclap/CmdLine.h"


int main (const int argc, const char*argv[]){
    std::string fasta_file_path;

    try {
        TCLAP::CmdLine cmd("", ' ', "0.1");

        // Read input fasta filename
        TCLAP::ValueArg<std::string> fasta_file_pathArg("F", "fasta_file_path",
                                                         "Input probe sequences (unpaired) files ", true, "",
                                                         "file1.fastq,file2.fastq", cmd);
        cmd.parse(argc, argv);

        fasta_file_path = fasta_file_pathArg.getValue();

    } catch (TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        return 1;
    }

    std::cout << fasta_file_path << std::endl;
    MarkerData md(fasta_file_path);

    std::cout << "First sequence: " << md.bases[0].ToString()<< std::endl;
    std::cout << md.probe_ids[0] << std::endl;

        return 0;
}