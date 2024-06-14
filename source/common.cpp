//
//  common.cpp
//
//  Copyright (c) 2013 Bruno Nevado. GNU license.
//

#include <fstream>
#include <algorithm>

#include "common.h"

// from IUPAC

vector <char> fromIUPAC(char in){
    vector <char> toreturn;
    
    if( in == 'a' || in == 'c' || in == 'g' || in == 't' || in == 'n' ){
        toreturn.push_back(in);
        toreturn.push_back(in);
        return toreturn;
    }
    else if ( in == 'm' ){
        toreturn.push_back('a');
        toreturn.push_back('c');
        return toreturn;
    }
    else if ( in == 'r' ){
        toreturn.push_back('a');
        toreturn.push_back('g');
        return toreturn;
    }
    else if ( in == 'w' ){
        toreturn.push_back('a');
        toreturn.push_back('t');
        return toreturn;
    }
    else if ( in == 's' ){
        toreturn.push_back('c');
        toreturn.push_back('g');
        return toreturn;
    }
    else if ( in == 'k' ){
        toreturn.push_back('g');
        toreturn.push_back('t');
        return toreturn;
    }
    else if ( in == 'y' ){
        toreturn.push_back('c');
        toreturn.push_back('t');
        return toreturn;
    }
    else{
        cerr << "ERROR (fromIUPAC): invalid code: " << in << "\n";
        exit(1);
    }
    
}

int sum_vec_ints ( vector <int> in ) {
    int sum = 0;
    for (unsigned int i = 0; i < in.size(); i++) {
        sum += in.at(i);
    }
    return(sum);
}


char toIUPAC(string in){
    //  cout << " called with " << in.at(0) << "\n";
    if( in.at(0) == in.at(1) ){
        return in.at(0);
    }
    else if ( (in.at(0) == 'a' && in.at(1) == 'c') || (in.at(1) == 'a' && in.at(0) == 'c') ){
        return 'm';
    }
    else if ( (in.at(0) == 'a' && in.at(1) == 'g') || (in.at(1) == 'a' && in.at(0) == 'g') ){
        return 'r';
    }
    else if ( (in.at(0) == 'a' && in.at(1) == 't') || (in.at(1) == 'a' && in.at(0) == 't') ){
        return 'w';
    }
    else if ( (in.at(0) == 'c' && in.at(1) == 'g') || (in.at(1) == 'c' && in.at(0) == 'g') ){
        return 's';
    }
    else if ( (in.at(0) == 't' && in.at(1) == 'g') || (in.at(1) == 't' && in.at(0) == 'g') ){
        return 'k';
    }
    else if ( (in.at(0) == 'c' && in.at(1) == 't') || (in.at(1) == 'c' && in.at(0) == 't') ){
        return 'y';
    }
    else{
        cerr << "ERROR (turnIUPAC): no code available for " << in.at(0) << in.at(1) << "\n";
        exit(1);
    }
    
    
}

unsigned int count_ignore_case ( char base, vector <char> DNAbases , vector <int> counts ){
    
    unsigned int count = 0;
    
    for (unsigned int i = 0; i < DNAbases.size() ; i++) {
        if ( DNAbases.at(i) == base || DNAbases.at(i) == toupper( base ) ) {
            count += counts.at(i);
        }
    }
    return count;
}

vector <char> get_uniques_lc ( vector <char> bases_to_check ){
    // from http://www.cplusplus.com/forum/general/14268/
    vector<char> result;
    vector<bool> duplicated(bases_to_check.size(), false);
    
    for (unsigned int i = 0; i < bases_to_check.size(); i++) {
        for (unsigned int j = i + 1 ; j < bases_to_check.size(); j++) {
            if (bases_to_check.at(i) == bases_to_check.at(j)) {
                duplicated.at(j) = true;
            }
        }
    }
    
    for (unsigned int i = 0; i < bases_to_check.size(); i++) {
        if (!duplicated.at(i)) {
            result.push_back(tolower(bases_to_check.at(i)));
        }
    }
    return result;
}

void msplit(string s , string delim,  vector<string> * output){
    
    //
    if( output->size() != 0 ){
        std::cerr << "WARNING: result vector for msplit is not empty!" << std::endl;
        for (unsigned int i = 0; i < output->size(); i++) {
            std::cerr << output->at(i) << " " ;
        }
        std::cerr << std::endl;
    }
    

    
    
    //
    unsigned long start = 0U;
    unsigned long end = s.find(delim);
    while (end != std::string::npos)
    {
        output->push_back(s.substr(start, end - start));
        start = end + delim.length();
        end = s.find(delim, start);
    }
    output->push_back(s.substr(start, end - start));
    
}

void msplit2( const string& s , string delim,  vector<string> * output){
    
    //
    if( output->size() != 0 ){
        std::cerr << "WARNING: result vector for msplit is not empty!" << std::endl;
        for (unsigned int i = 0; i < output->size(); i++) {
            std::cerr << output->at(i) << " " ;
        }
        std::cerr << std::endl;
    }
    

    //
    unsigned long start = 0U;
    unsigned long end = s.find(delim);
    while (end != std::string::npos)
    {
        output->push_back(s.substr(start, end - start));
        start = end + delim.length();
        end = s.find(delim, start);
    }
    output->push_back(s.substr(start, end - start));
    
}


vector <float> calc_dist( vector <int> in ){
    // vector <float> toreturn(2,0);  //mean var
    float mean = 0, var = 0;
    int count = 0;
    // starts at 1 to exclude 0 class (whre there is no missing data)
    for (unsigned int i = 1; i < in.size(); i++) {
        mean += ( in.at(i)*i );
        count += in.at(i);
    }
    mean = mean/count;
    
    for (unsigned int i = 1; i < in.size(); i++) {
        var += (i - mean) * (i - mean) * in.at(i);
        
    }
    
    var = var / (count - 1);
    vector <float> toreturn;
    toreturn.push_back(mean);
    toreturn.push_back(var);
    return toreturn;
}


void check_regions ( string regions, vector < unsigned int > &starts, vector < unsigned int > &ends, vector <string> &names ) {
    
    // case 1: single region defined in command line
    if( ( count(regions.begin(), regions.end(), ':')) ){
        // single region
        vector <string> fields;
        msplit(regions, ":", &fields );
        starts.push_back(atoi(fields.at(0).c_str()));
        ends.push_back(atoi(fields.at(1).c_str()));
        names.push_back("NA");
        if( starts.at(0) == 0 ||  ends.at(0) == 0){
            cout << "ERROR (regions) : error in definition of region to output (received " << regions << ")\n";
            exit(1);
        }
        
        
    }
    else{
        string line;
        ifstream infile;
        
        infile.open(regions.c_str());
        
        
        if( !infile.is_open() ) {
            cout << "ERROR (regions): Unable to open for output regions file " << regions << endl;
            exit(1);
        }
        while (getline(infile,line)) {
            if (line.length() == 0) {
                continue;
            }
            
            vector <string> fields;
            msplit( line, "\t", &fields );
            
            if ( fields.size() != 3 && fields.size() != 2 ){
                
                cout << "ERROR (regions): Malformatted region in file " << regions << ": " << line <<  endl;
                exit(1);
                
            }
            
            if( fields.size() == 3 ){
                starts.push_back(atoi(fields.at(0).c_str()));
                ends.push_back(atoi(fields.at(1).c_str()));
                names.push_back(fields.at(2));
                
            }else{
                starts.push_back(atoi(fields.at(0).c_str()));
                ends.push_back(atoi(fields.at(1).c_str()));
                names.push_back("NA");
                
            }
            
            
            
            
        }
        
        infile.close();
        if (starts.size() == 0 ) {
            cout << "ERROR (regions): no regions defined in input file " << regions << "?\n";
            exit(1);
        }
        
        
        for (unsigned int k = 0; k < starts.size(); k++) {
            if( starts.at(k) == 0 ||  ends.at(k) == 0 || ends.at(k) <= starts.at(k) ){
                cout << "ERROR (regions) : error in definition of region number " << k + 1 << " in input file " << regions << endl;
                exit(1);
            }
        }
        
    }
    
}


bool contains_chars_not_in_string( std::string test, std::string allowed){
    std::size_t found;
    for (unsigned int i = 0; i < test.size(); i++) {
        found=allowed.find(test.at(i));
        if (found==std::string::npos)
            return true;
        
    }
    return false;
}


bool is_abba(std::string in ){
    
    if( in.at(0) == in.at(3) && in.at(1) == in.at(2) && in.at(0) != in.at(1) ){
        return true;
    }
    else{
        return false;
    }
    
}


bool is_baba(std::string in){
    
    if( in.at(0) == in.at(2) && in.at(1) == in.at(3) && in.at(0) != in.at(1) ){
        return true;
    }
    else{
        return false;
    }
    
}


