//
//  main.cpp
//  fasta2paml
//
//  Created by brunonevado on 10/10/2014.
//  Copyright (c) 2014 ___BRUNO_DOPS___. All rights reserved.
//


#include <iostream>
#include <fstream>
#include <sstream>

#include "fasta.h"
#include "common.h"
#include "args.h"

#include <Bpp/Phyl.all> /* this includes all models */


// 20.03.2015: added different cut off values for sites and sequences
// 23.04.2015: added strictNames option
// 12.05.2015: added compatibility for clade names in intree

void help(){
    std::cout << "###################\n  fasta2paml2 12052015 \n###################" << std::endl;;
    std::cout << "Makes phylip file from fasta file after:" << std::endl;
    std::cout << "      - removing last codon if it is a stop" << std::endl;
    std::cout << "      - masking any stop codons within alignment" << std::endl;
    std::cout << "      - removing sequences with ambiguous bases % above threshold (-cutSeqs)" << std::endl;
    std::cout << "      - removing codon positions with ambiguous proportion above threshold (-cutCodons)" << std::endl;
    
    std::cout << "Usage: fasta2paml -infile in.fas -outfile out_prefix -names names.txt -tree tree.tre -cut n -dropBranches 1/0 -unroot 1/0 -strictNames 1/0" << std::endl;
    std::cout << "Args: -infile: fasta file with sequences." << std::endl;
    std::cout << "      -outfile: prefix for output files (will write to outfile.phy and outfile.tre)." << std::endl;
    std::cout << "      -names: name of sequences to write - sequences or leaves not named will be removed." << std::endl;
    std::cout << "      -tree: input tree file, newick format, with or without branch lengths." << std::endl;
    std::cout << "      -cutSeqs|cutCodons: n threshold to remove sequences/codons (e.g. with 50, codons where more than 50% individuals have missing data, or sequences with > 50% missing bases, are removed)." << std::endl;
    std::cout << "      -dropBranches: whether to keep (0) or remove (1) branch lengths if present in tree file." << std::endl;
    std::cout << "      -unroot: whether to unroot the tree before printing (otherwise the tree is printed as read, after removing tips)." << std::endl;
    std::cout << "      -strictNames: if set to 1, species in '-names' missing in alignments will raise error, otherwise (set to 0) will raise warning only." << std::endl;

    
    std::cout << "Only sequences listed in file names.txt are printed (1 sequence name per line)." << std::endl;
    std::cout << "Using '-names all' will print all sequences in infile." << std::endl;
    std::cout << "Sequences removed due to missing data are also removed from tree file." << std::endl;
    std::cout << "Sequences in tree that are not present in names file (or fasta file) are also removed from tree." << std::endl;

    
};

int main(int argc, const char * argv[])
{
    
     //get args
    ////////////////////
    sargs myargs;
    try{
        myargs = args::getargs(argc, argv, std::vector<std::string> {"infile","outfile","names","tree"}, std::vector<std::string> {"dropBranches","unroot","strictNames"}, std::vector<std::string>  {"cutSeqs","cutCodons"}, std::string {}, std::string {}); }
    catch (std::string e){
        std::cout << " Args failed: " << e << std::endl;
        help();
        exit(1);
    }
    
    
    
    std::string infile = myargs.args_string.at(0);
    std::string outfile = myargs.args_string.at(1);
    std::string names = myargs.args_string.at(2);
    std::string treefile = myargs.args_string.at(3);
    
    int maxCutSeqs = myargs.args_int.at(0);
    int maxCutCodons = myargs.args_int.at(1);
    bool dropBlens = myargs.args_booleans.at(0);
    bool unroot = myargs.args_booleans.at(1);
    bool strict_names = myargs.args_booleans.at(2);

    
    int nStopsmasked = 0;
    int nSeqsRemoved = 0;
    int nCodonsRemoved = 0;
    
    
    
    // READ TREE FILE
    bpp::Newick reader;

    reader.enableExtendedBootstrapProperty("cladeNames");
    
    bpp::Tree  * atree = reader.read(treefile);
    bpp::TreeTemplate<bpp::Node> tree1 (*atree);
    
    
    // READ FASTA FILE
    fasta afasta(1);
    afasta.read_fasta_file(infile);
    fasta workFasta(1);
    
    // subset fasta if -names
    if( names != "all" ){
        // extract subset of sequences
        std::vector < std::string> seq_names;
        ifstream fh_names(names);
        if( !fh_names.is_open() ){
            std::cerr << "<fasta2paml> ERROR: unable to open for reading names file " << names << std::endl;
            exit(1);
        }
        std::string cline;
        while (getline(fh_names, cline )) {
            seq_names.push_back(cline);
        }
        //fasta newfasta(afasta.num_lines());
        workFasta.new_fasta_from_names(afasta, seq_names, strict_names);
    }
    else{
        workFasta = afasta;
    }
    
    
    // here check tree against names
    std::vector <std::string > namesInTree = tree1.getLeavesNames();
    std::vector <std::string > namesInFasta = workFasta.show_names();
    for (unsigned int ileave = 0; ileave < namesInTree.size() ; ileave++) {
        if(std::find(namesInFasta.begin(), namesInFasta.end(), namesInTree.at(ileave)) == namesInFasta.end()){
            //std::cout << "removing leave " << ileave << " is " << namesInTree.at(ileave) << std::endl;
            bpp::TreeTemplateTools::dropLeaf(tree1, namesInTree.at(ileave) );
        }
        

    }
    
    // drop branch lengths if requested
    if( dropBlens ){
        for(auto inode : tree1.getNodes()){
            inode->deleteDistanceToFather();
        }
     }
    
    
   // exit(0);
    // remove and mask stop codons
    workFasta.remove_last_if_stop();
    nStopsmasked = workFasta.mask_stop_codons();
    if( nStopsmasked != 0 ){
        std::cout << "<fasta2paml> WARNING: "<< nStopsmasked << " STOP CODON(S) FOUND AND MASKED IN MIDDLE OF ALIGNMENT IN FILE " << infile << std::endl;
    }
    
    // remove sequences if too many ambiguous sites (also trim tree)
    std::vector <int> seqs_to_remove;
    for (unsigned int iseq = 0; iseq < workFasta.num_lines(); iseq++) {
        float missing = 100.0 * float(workFasta.get_num_amb(iseq))/float(workFasta.num_bases());
        if(  missing > maxCutSeqs ){
            seqs_to_remove.push_back(iseq);
        }
    }
    nSeqsRemoved = int (seqs_to_remove.size());
    if( nSeqsRemoved != 0 ){
        std::cout << "<fasta2paml> WARNING: "<< nSeqsRemoved << " SEQUENCES REMOVED DUE TO MISSING DATA IN FILE " << infile << std::endl;
    }
    if( nSeqsRemoved == workFasta.num_lines() ){
        std::cerr << "<fasta2paml> ERROR: ALL SEQUENCES HAVE MISSING PROPORTION ABOVE THRESHOLD!"  << std::endl;
        exit(1);
    }
    for( int i = seqs_to_remove.size() -1; i >= 0; i-- ){
        std::cout << "dropping " << seqs_to_remove.at(i) << " which is " << workFasta.show_name(seqs_to_remove.at(i)) << std::endl;
        bpp::TreeTemplateTools::dropLeaf(tree1, workFasta.show_name(seqs_to_remove.at(i)) );
        workFasta.delete_seq(seqs_to_remove.at(i));
    }
    
    // now codons
    nCodonsRemoved = workFasta.remove_ambiguous_codons(maxCutCodons);
    
    std::clog << "<fasta2paml> SUMMARY (stopsMasked, seqsRemoved, codonsRemoved)\t" << infile << "\t" << nStopsmasked << "\t" << nSeqsRemoved <<  "\t" << nCodonsRemoved << std::endl;
    
    std::stringstream outPhy;
    outPhy << outfile << ".paml.phy";
    workFasta.write_phylip( outPhy.str() );
    
    if(unroot && tree1.isRooted()){
        tree1.unroot();
    }
    
    std::stringstream outTre;
    outTre << outfile << ".paml.tre";
    reader.write( tree1 , outTre.str());

    return 0;
}

