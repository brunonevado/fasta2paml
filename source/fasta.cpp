//
//  fasta.cpp
// Created by Bruno Nevado on 17/06/2013.
//  Copyright (c) 2013 Bruno Nevado. All rights reserved.
//
#include "fasta.h"


#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include <locale>

#include "fasta.h"

// FASTA OBJECT CONSTRUCTOR

fasta::fasta ( int num_inds, int len ){
    if(len == 0){
        names.reserve(num_inds);
        matrix.reserve(num_inds);
    }
    else{
        matrix.resize(num_inds);
        names.resize(num_inds);
        
        for(int i = 0; i < num_inds; i++)
            matrix.at(i).resize(len);
        /*
         string newline (len , NULL);
         for(int i = 0; i < num_inds; i++)
         matrix.push_back(newline);
         
         
         */
    }
    
        //acgtRY SWKM
        //  tgcayr swmk
        //genetic_code. wikipedia... oh boy ... http://en.wikipedia.org/wiki/DNA_codon_table
        genetic_code.insert(make_pair("gct","A"));
        genetic_code.insert(make_pair("gcc","A"));
        genetic_code.insert(make_pair("gca","A"));
        genetic_code.insert(make_pair("gcg","A"));
        
        genetic_code.insert(make_pair("cgt","R"));
        genetic_code.insert(make_pair("cgc","R"));
        genetic_code.insert(make_pair("cga","R"));
        genetic_code.insert(make_pair("cgg","R"));
        genetic_code.insert(make_pair("aga","R"));
        genetic_code.insert(make_pair("agg","R"));
        
        genetic_code.insert(make_pair("aat","N"));
        genetic_code.insert(make_pair("aac","N"));
        
        genetic_code.insert(make_pair("gat","D"));
        genetic_code.insert(make_pair("gac","D"));
        
        genetic_code.insert(make_pair("tgt","C"));
        genetic_code.insert(make_pair("tgc","C"));
        
        genetic_code.insert(make_pair("caa","Q"));
        genetic_code.insert(make_pair("cag","Q"));
        
        genetic_code.insert(make_pair("gaa","E"));
        genetic_code.insert(make_pair("gag","E"));
        
        genetic_code.insert(make_pair("ggt","G"));
        genetic_code.insert(make_pair("ggc","G"));
        genetic_code.insert(make_pair("gga","G"));
        genetic_code.insert(make_pair("ggg","G"));
        
        genetic_code.insert(make_pair("cat","H"));
        genetic_code.insert(make_pair("cac","H"));
        
        genetic_code.insert(make_pair("att","I"));
        genetic_code.insert(make_pair("atc","I"));
        genetic_code.insert(make_pair("ata","I"));
        
        genetic_code.insert(make_pair("tta","L"));
        genetic_code.insert(make_pair("ttg","L"));
        genetic_code.insert(make_pair("ctt","L"));
        genetic_code.insert(make_pair("ctc","L"));
        genetic_code.insert(make_pair("cta","L"));
        genetic_code.insert(make_pair("ctg","L"));
        
        genetic_code.insert(make_pair("aaa","K"));
        genetic_code.insert(make_pair("aag","K"));
        
        genetic_code.insert(make_pair("atg","START"));
        
        genetic_code.insert(make_pair("ttt","F"));
        genetic_code.insert(make_pair("ttc","F"));
        
        genetic_code.insert(make_pair("cct","P"));
        genetic_code.insert(make_pair("ccc","P"));
        genetic_code.insert(make_pair("cca","P"));
        genetic_code.insert(make_pair("ccg","P"));
        
        genetic_code.insert(make_pair("tct","S"));
        genetic_code.insert(make_pair("tcc","S"));
        genetic_code.insert(make_pair("tca","S"));
        genetic_code.insert(make_pair("tcg","S"));
        genetic_code.insert(make_pair("agt","S"));
        genetic_code.insert(make_pair("agc","S"));
        
        genetic_code.insert(make_pair("act","T"));
        genetic_code.insert(make_pair("acc","T"));
        genetic_code.insert(make_pair("aca","T"));
        genetic_code.insert(make_pair("acg","T"));
        
        genetic_code.insert(make_pair("tgg","W"));
        
        genetic_code.insert(make_pair("tat","Y"));
        genetic_code.insert(make_pair("tac","Y"));
        
        genetic_code.insert(make_pair("gtt","V"));
        genetic_code.insert(make_pair("gtc","V"));
        genetic_code.insert(make_pair("gta","V"));
        genetic_code.insert(make_pair("gtg","V"));
        
        genetic_code.insert(make_pair("taa","STOP"));
        genetic_code.insert(make_pair("tga","STOP"));
        genetic_code.insert(make_pair("tag","STOP"));
        genetic_code.insert(make_pair("tar","STOP"));
        genetic_code.insert(make_pair("tra","STOP"));


}

// FASTA OBJECT CALCULATE ALLELE COUNT

int fasta::calc_dac ( char alt , unsigned int site ) {
    int freq = 0;
    
    for(int ind = 1; ind < int ( matrix.size() ); ind++){
        if( matrix[ind].at(site-1) == alt){
            freq++;
        }
    }
    return freq;
}

// FASTA OBJECT OUTPUT SNP FILES


void fasta::outputSNPfiles( int cind, string prefix ){
    
    ofstream streamHomo, streamHetero;
    
    stringstream cindstr;
    cindstr << cind;
    string outfileHomo = prefix + ".Ind" + cindstr.str() + ".fixed";
    string outfileHetero = prefix + ".Ind" + cindstr.str() + ".snps";
    
    streamHomo.open(outfileHomo.c_str());
    if( !streamHomo.is_open() ){
        cerr << "ERROR (outputSNPfiles): unable to open for output file " << outfileHomo << "\n";
        exit(1);
    }
    streamHetero.open(outfileHetero.c_str());
    if( !streamHetero.is_open() ){
        cerr << "ERROR (outputSNPfiles): unable to open for output file " << outfileHetero << "\n";
        exit(1);
        
    }
    int i = (cind -1) * 2 + 1;
    for( unsigned int j = 0; j < this->num_bases(); j++ ){
        
        if(matrix[0][j] == matrix[i][j] && matrix[i][j] == matrix[i+1][j])
            continue;
        // Fixed diff
        if ( matrix[0][j] != matrix[i][j] && matrix[i][j] == matrix[i+1][j] ){
            streamHomo << j+1 << "\t" << matrix[0][j] << "\t" << matrix[i][j] << "\t" << matrix[i+1][j]
            << "\t" << calc_dac( matrix[i][j] , j+1 ) << "\n";
            
        }
        // heterozygous SNP
        else if ( (matrix[0][j] != matrix[i][j] && matrix[0][j] == matrix[i+1][j] )
                 || (matrix[0][j] == matrix[i][j] && matrix[0][j] != matrix[i+1][j]) ){
            streamHetero << j+1 << "\t" << matrix[0][j] << "\t" << matrix[i][j] << "\t" << matrix[i+1][j]
            << "\t" << calc_dac( ( matrix[0][j] == matrix[i+1][j] ?  matrix[i][j] :  matrix[i+1][j] ) , j+1 ) << "\n";
        }
        else{
            cout << "WARNING: position " << j+1 << " is triallelic (ref:" << matrix[0][j] << ",individual " << cind << ":" << matrix[i][j] << matrix[i+1][j] << ")\n" ;
        }
        
    }
    streamHomo.close();
    streamHetero.close();
    
    
}


// FASTA CLASS SPEAKING FUNCTION

void fasta::info_to_stdout () {
    
    clog << "Infile: " << this->input_file() << " (" << this->num_lines() << " sequences, " <<
    this->num_bases() << " bp long)\n";
    
    /*
     cout << "Name of sequences: ";
     for(int i = 0; i < this->num_lines(); i++)
     cout << this->names[i] << " ";
     cout << "\n";
     
     cout<< "And the data matrix is:\n";
     for(int i = 0; i < matrix.size(); i++)
     cout << matrix[i] << "\n";
     
     
     */
}

// FASTA CLASS FUNCTION: READ FASTA FILE

void fasta::read_fasta_file( std::string infas){
    int cind = 0;
    matrix.clear();
    std::locale loc;
    std::string line;
    std::ifstream infile_fas (infas.c_str());
    if (infile_fas.is_open())
    {
        infile = infas;
        while ( ! infile_fas.eof() )
        {
            getline(infile_fas, line);
            
            if( line[0] == '>' ){
                cind++;
                matrix.resize(cind);
                std::string name = line.substr(1);
                names.push_back(name);
            }
            else {
                matrix.at(cind-1).append( line );
            }
        }
    }
    
    else {
        std::cerr << "ERROR (read_fasta_file): Unable to open infile " << infas << "\n" ;
        exit(1);
    }
    infile_fas.close();
    
    for (unsigned int l = 0; l < matrix.size(); l++) {
        for (unsigned s = 0; s < matrix.at(l).length(); s++) {
            matrix.at(l).at(s) = tolower(matrix.at(l).at(s), loc);
        }
    }
    
}

// CHECK WHETHER ALL SEQS HAVE SAME LENGTH
int fasta::is_aligned () {
    unsigned int previous_len = 0;
    for(unsigned int i = 0; i < matrix.size() ; i++){
        unsigned int current_len = 0;
        for (unsigned int j = 0; j <  matrix.at(i).size() ; j++) {
            if ( matrix[i][j] != char (NULL)  ){
                current_len++;
            }
            
        }
        if( i == 0){
            previous_len = current_len;
            continue;
        }
        else if (current_len != previous_len){
            return (i);
        }
    }
    
    return 0;
}

// WRITES RANDOM SEQUENCE TO FILE

void fasta::write_random( string out, string name, unsigned int len ){
    
    
    ofstream outputFile;
    char bases[4] = { 'A', 'C', 'G', 'T' };
    
    outputFile.open(out.c_str());
    
    if( !outputFile.is_open() ){
        cerr << "ERROR (write_random): unable to open for output file " << out << "\n";
        exit(1);
    }
    
    outputFile << ">" << name << "\n";
    
    for( unsigned int site = 0; site < len; site++){
        outputFile << bases[ rand() % 4 ];
    }
    outputFile << "\n";
    outputFile.close();
    
}

// WRITE ALIGNMENT TO FILE

void fasta::write_to_file( string out, int append ){
    locale loc;
    ofstream outputFile;
    
    if( append == 0 )
        outputFile.open(out.c_str());
    else
        outputFile.open(out.c_str(), ios::app );
    
    if( !outputFile.is_open() ){
        cerr << "ERROR (write_to_file): unable to open for output file " << out << "\n";
        exit(1);
    }
    if( names.size() == matrix.size() ){
        for( unsigned int i = 0; i < matrix.size(); i++){
            outputFile << ">" << names[i] << "\n";
            for (unsigned int site = 0; site < matrix.at(i).length(); site++) {
                outputFile << toupper(matrix[i][site],loc); //  matrix[i][site];
            }
            outputFile << endl;
        }
    }
    else{
        for( unsigned int i = 0; i < matrix.size(); i++){
            outputFile << ">i" << i+1 << "\n"; //
            for (unsigned int site = 0; site < matrix.at(i).length(); site++) {
                outputFile << toupper(matrix[i][site],loc); //  matrix[i][site];
            }
            outputFile << endl;
        }
    }
    outputFile.close();
    
}

void fasta::remove_multi_hits(){
    
    vector <unsigned int> mh_pos;
    
    for ( unsigned int site = 0; site < this->num_bases(); site++){
        char ref = this->matrix[0][site] , alt = 'N';
        for( unsigned int ind = 1; ind < this->num_lines(); ind++ ){
            if( matrix[ind][site] != ref ){
                if (alt == 'N'){
                    alt = matrix[ind][site];
                }
                else if ( matrix[ind][site] != alt ){
                    mh_pos.push_back(site);
                    
                }
            }
        }
        
    }
    
    if( mh_pos.size() > 0 ){
        cout << "  Multiple hits found: " << mh_pos.size()
        << ". Sites masked as invariable: ";
        for( unsigned long mh = 0; mh < mh_pos.size(); mh++ ){
            cout << " " << mh_pos[mh]+1;
            for( unsigned int ind = 1; ind < matrix.size(); ind++){
                matrix[ind][mh_pos[mh]] = matrix[0][mh_pos[mh]];
                
            }
        }
        cout << "\n";
    }
    else{
        cout << "  No multiple hits found\n";
    }
    
}


// FOR KICKS

void fasta::forkicks(){
    ofstream outputFile;
    outputFile.open("table.txt");
    
    double sites=0;
    double start = 1;
    for( unsigned int site = 0; site < this->num_bases(); site++){
        int present=0;
        for( unsigned int seq = 0; seq < this->num_lines(); seq++){
            
            if(matrix[seq][site] != 'n'){
                present++;
            }
            if(present == int ( this->num_lines() ) )
                sites++;
            
        }
        
        if(( site % 500000) == 0 ){
            outputFile << fixed << start << "\t" << site << "\t" << sites << "\n";
            start=site;
            sites=0;
            
        }
        
        
    }
    
    outputFile.close();
    
}


int fasta::mask_missing_in_last(){
    int nsites = 0;
    
    for( unsigned int isite = 0; isite < matrix.at(0).size(); isite++ ){
        if( matrix.at(matrix.size()-1).at(isite) == 'n' || matrix.at(matrix.size()-1).at(isite) == 'N' ){
            nsites++;
            for( int iline = 0; iline < matrix.size() -1 ; iline++){
                matrix.at(iline).at(isite) = 'n';
            }
        }
    }
    return nsites;
}


// FASTA ADD SEQUENCE OPTION

void fasta::append_seq( fasta in, bool fill ){
    
    for( unsigned int iline = 0; iline < matrix.size(); iline++ ){
        
        for( unsigned int isite = matrix.at(iline).length(); isite < in.num_bases(); isite++ ){
            matrix.at(iline).append("n");
        }
    }
    matrix.push_back(in.matrix.at(0));
    names.push_back( in.names.at(0));
    
}

void fasta::resize_matrix( unsigned int start, unsigned int end ){
    for( unsigned int line = 0; line < matrix.size(); line++ ){
        string sub = matrix.at(line).substr(start - 1, end - start + 1  );
        matrix.at(line) = sub;
    }
    
}

void fasta::mask_base ( unsigned int site1, char newchar, int index1){
    
    if ( site1 == 0 ) {
        cerr << "ERROR (fasta::mask_base): site number should be 1-based" << endl;
        exit(1);
    }
    if( index1 == -1 ){
        for (unsigned int line = 0; line < matrix.size(); line++) {
            
            if (site1  > matrix.at(line).size() ) {
                cerr << "ERROR (fasta::mask_base): site to mask (" <<  site1 <<") out of bounds" << endl;
                exit(1);
            }
            else{
                matrix.at(line).at( site1 - 1 ) = tolower ( newchar );
            }
        }
        
    }
    else{
        if (site1  > matrix.at(index1 -1 ).size() ) {
            cerr << "ERROR (fasta::mask_base): site to mask (" <<  site1 <<") out of bounds" << endl;
            exit(1);
        }
        else{
            matrix.at(index1 -1).at( site1 - 1 ) = tolower ( newchar );
        }
    }
}


void fasta::new_fasta_from_inds ( fasta &infasta, vector <int> index0 ){
    
    infile = infasta.infile;
    for (unsigned int line = 0; line < index0.size() ; line++ ) {
        matrix.push_back( infasta.matrix.at(index0.at(line)));
        names.push_back( infasta.names.at(index0.at(line)));
    }
}


void fasta::new_fasta_from_names(fasta &infasta, vector<std::string> inames, bool strict){
    
    infile = infasta.infile;
    
    
    
    for (unsigned int iname = 0; iname < inames.size() ; iname++ ) {
        bool found = false;
        for (unsigned int iline = 0; iline < infasta.num_lines() ; iline++ ) {
            
            if( infasta.names.at(iline) == inames.at(iname) ){
                matrix.push_back(infasta.matrix.at(iline));
                names.push_back(infasta.names.at(iline));
                found = true;
                break;
            }
        }
        if(!found and strict){
            std::cerr << "ERROR in fasta_from_names: sequence " << inames.at(iname) << " not found in file " << infasta.infile << std::endl;
            exit(1);
        }
        else if (!found){
            std::cerr << "WARNING in fasta_from_names: sequence " << inames.at(iname) << " not found in file " << infasta.infile << std::endl;

        }
    }
}




int fasta::concatenate_alignments( fasta & al2){
    if( matrix.size() != al2.num_lines()){
        cerr << "ERROR (concatenate fasta): number of sequences do not match!"  << endl;
        exit(1);
    }
    
    for ( unsigned int i = 0; i < matrix.size(); i++ ) {
        matrix.at(i).append( al2.show_seq(i) );
    }
    return 0;
}


void fasta::compare_species ( int wlen ){
    
    int bases_seen = 0, bases_differ = 0, bases_in1 = 0, bases_in2 = 0;
    int nwindows = (wlen != 0) ? this->num_bases()/wlen : 1;
    if (wlen != 0 && this->num_bases() % wlen != 0) {
        nwindows++;
    }
    std::cout << "Window\tStart\tEnd\tSP1_bases\tSP2_bases\tOverlap\tP_DIFF" << std::endl;
    
    
    for (int n = 0; n < nwindows; n++) {
        int start = n*wlen+1;
        unsigned int end = (n*wlen+ wlen <= this->num_bases() && wlen != 0 )? n*wlen+ wlen : this->num_bases() ;
        
        for (unsigned int site = start; site <= end; site++) {
            if( matrix.at(0).at(site -1 ) == 'a'
               || matrix.at(0).at(site -1 ) == 'c'
               || matrix.at(0).at(site -1 ) == 'g'
               || matrix.at(0).at(site -1 ) == 't'
               ){
                bases_in1++;
            }
            if( matrix.at(1).at(site -1 ) == 'a'
               || matrix.at(1).at(site -1 ) == 'c'
               || matrix.at(1).at(site -1 ) == 'g'
               || matrix.at(1).at(site -1 ) == 't'
               ){
                bases_in2++;
            }
            
            if(( matrix.at(0).at(site -1 ) == 'a'
                || matrix.at(0).at(site -1 ) == 'c'
                || matrix.at(0).at(site -1 ) == 'g'
                || matrix.at(0).at(site -1 ) == 't'
                )
               &&
               ( matrix.at(1).at(site -1 ) == 'a'
                || matrix.at(1).at(site -1 ) == 'c'
                || matrix.at(1).at(site -1 ) == 'g'
                || matrix.at(1).at(site -1 ) == 't'
                )){
                   bases_seen++;
                   
                   if(matrix.at(0).at(site -1 ) != matrix.at(1).at(site -1 )){
                       bases_differ++;
                   }
                   
                   
               }
        }
        
        std::cout << n+1 << "\t" << start << "\t" << end << "\t" << bases_in1 << "\t" << bases_in2  << "\t" <<bases_seen << "\t" << float(bases_differ)/float(bases_seen) << std::endl;
        bases_seen = 0;
        bases_differ = 0;
        bases_in1 = 0;
        bases_in2 = 0;
        
    }
    
    
}

void fasta::build_from_fasta(fasta & in, std::vector <int> sites){
    matrix.resize(in.num_lines());
    names = in.names;
    for(unsigned int isite = 0; isite < sites.size(); isite++ ){
        
        int site_to_push =sites.at(isite) - 1;
        for( unsigned int iline = 0; iline < matrix.size(); iline++ ){
            matrix.at(iline).push_back( in.matrix.at(iline).at(site_to_push) );
        }
        
    }
}



void fasta::write_phylip(std::string outfile, bool append){
    
    std::ofstream fh_out;
    
    if( append  )
        fh_out.open(outfile.c_str(), ios::app);
    else
        fh_out.open(outfile.c_str() );
    
    
    if( !fh_out.is_open() ){
        cerr << "ERROR (fasta2phylip): unable to open for output file " << outfile << std::endl;
        exit(1);
    }
    fh_out << this->num_lines() << " " << this->num_bases() << std::endl;
    for (unsigned int iind = 0; iind < this->num_lines(); iind++) {
        std::string iname = names.at(iind).substr(0,10);
        iname.insert(iname.end(), 10-iname.length(), ' ' );
        fh_out << iname << " " << matrix.at(iind) << std::endl;
    }
    fh_out.close();
}


void fasta::clean_matrix(){
    for ( int isite = this->num_bases()-1 ; isite >= 0 ; isite-- ) {
        std::string column;
        for (unsigned int iline = 0; iline < this->num_lines(); iline++) {
            column.push_back(matrix.at(iline).at(isite));
        }
        
        if( contains_chars_not_in_string(column, "acgt") ){
            for (unsigned int iline = 0; iline < this->num_lines(); iline++) {
                matrix.at(iline).erase( matrix.at(iline).begin()+isite );
            }
        }
        
        
    }

}

int fasta::get_num_vars(){
    int nvars = 0;
    for ( int isite = this->num_bases()-1 ; isite >= 0 ; isite-- ) {
        std::string column, pattern;
        for (unsigned int iline = 1; iline < this->num_lines(); iline++) {
            column.push_back(matrix.at(iline).at(isite));
        }
        pattern.push_back(matrix.at(0).at(isite));
        
        if( contains_chars_not_in_string(column, pattern )){
            nvars++;
        }
        
        
    }
    return nvars;
}

int fasta::iupac_to_haploid_random(){
    
    int counter = 0;
    for (unsigned int iline = 0; iline < this->num_lines(); iline++ ) {
        for (unsigned int isite = 0; isite < this->matrix.at(iline).length(); isite++) {
            
            
            if( this->matrix.at(iline).at(isite) != 'a' && this->matrix.at(iline).at(isite) != 'c'
               && this->matrix.at(iline).at(isite) != 'g' && this->matrix.at(iline).at(isite) != 't'
               && this->matrix.at(iline).at(isite) != 'n' ){
                
                //std::cout <<this->matrix.at(iline).at(isite) << " => ";
                
                this->matrix.at(iline).at(isite) = fromIUPAC(this->matrix.at(iline).at(isite)).at(rand()%2);
                counter++;
                //std::cout <<this->matrix.at(iline).at(isite) << std::endl;
            }
            
            
        }
    }
    
    
    
    return counter;
    
}


void fasta::remove_last_if_stop(){
    bool remove = false;

    for (unsigned int i = 0; i < matrix.size(); i++) {
        
        std::string lastcodon = matrix.at(i).substr( matrix.at(i).length() - 3 );
        //std::cout << " INFILE " << infile << " last " << lastcodon
      //  << std::endl;

        if( genetic_code.count(lastcodon) == 0 ){ continue;}

        if( genetic_code.at(lastcodon) == "STOP"  ){
            remove = true;
            //std::cout << matrix.at(i) << std::endl;
            break;
        }

    }
    
    if(remove){
        for (unsigned i = 0; i < matrix.size(); i++) {
            matrix.at(i).erase(matrix.at(i).length()- 3 , matrix.at(i).length());
            //std::cout << matrix.at(i) << std::endl;
        }
    }


}


bool fasta::has_stop_codon(  ){
    
    for (unsigned int iline = 0; iline < matrix.size(); iline++) {
        for (unsigned int isite = 0; isite < matrix.at(iline).size(); isite += 3 ) {
            std::string codon = matrix.at(iline).substr( isite, 3 );
            
            if( genetic_code.count(codon) == 0 ){ continue;}
            //std::cout << codon << ":" << genetic_code.at(codon) << std::endl;
            
            
            if( genetic_code.at(codon) == "STOP"  ){return true;}
        }
    }
    
    return false;
    
}

int fasta::mask_stop_codons() {
    int masked = 0;
    for (unsigned int iline = 0; iline < matrix.size(); iline++) {
        for (unsigned int isite = 0; isite < matrix.at(iline).size(); isite += 3 ) {
            std::string codon = matrix.at(iline).substr( isite, 3 );
            
            if( genetic_code.count(codon) == 0 ){ continue;}
            //std::cout << codon << ":" << genetic_code.at(codon) << std::endl;
            
            
            if( genetic_code.at(codon) == "STOP"  ){
                matrix.at(iline).replace( isite, 3,"NNN");
                masked++;
                
            }
        }
    }
    
    return masked;
    
}

int fasta::get_num_amb( int ind0 ){
    int amb = 0;
    for (unsigned int isite = 0; isite < this->num_bases(); isite++ ) {
        if( matrix.at(ind0).at(isite) != 'a' &&   matrix.at(ind0).at(isite) != 'c'   &&   matrix.at(ind0).at(isite) != 'g' &&   matrix.at(ind0).at(isite) != 't' ){
            amb++;
        }
    }
    
    return amb;
}

void fasta::delete_seq(int ind0){
    matrix.erase(matrix.begin() + ind0);
    names.erase(names.begin() + ind0);
    
}


int fasta::remove_ambiguous_codons (int thr){
    int removed = 0;
    int siteCount = 0;
    int siteTot = this->num_bases();
    for (int isite = this->num_bases() - 3 ; isite >= 0 ; isite -= 3 ) {
        siteCount++;
        bpp::ApplicationTools::displayGauge(siteCount, siteTot,  '.', " Parsing codons " );

        int amb = 0;
        for (unsigned int iline = 0; iline < matrix.size(); iline++) {
            std::string codon = matrix.at(iline).substr( isite, 3 );
            if( genetic_code.count(codon) == 0 ){ amb++;}
            if ( 100.0 * float(amb)/float(matrix.size()) > thr ){
                break;
            }
        }
        //std::cout << this->num_bases() << " : " << float(amb)/float(matrix.size()) << ", thr: " << thr << std::endl ;

        if( 100.0 * float(amb)/float(matrix.size()) > thr ){
            removed++;
            for (unsigned int iline = 0; iline < matrix.size(); iline++) {
               matrix.at(iline).erase( matrix.at(iline).begin() + isite,  matrix.at(iline).begin() + isite + 3  );
            }
        }
        
    }
    return removed;
}
