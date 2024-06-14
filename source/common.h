//
//  common.h
//
//  Copyright (c) 2013 Bruno Nevado. GNU license.
//


#include <iostream>
#include <vector>
#include <string>
#include <algorithm>


using namespace std;

void msplit( string tosplit , string delim,  vector<string> * output );
void msplit2( const string& tosplit , string delim,  vector<string> * output );

// returns mean and var from a distribution
vector <float> calc_dist( vector <int> );

int sum_vec_ints ( vector <int> in ) ;

void check_regions ( string regions, vector <unsigned int> &starts, vector <unsigned int> &ends, vector <string> &names );

char  toIUPAC (string);
vector <char> fromIUPAC (char);


unsigned int count_ignore_case ( char base, vector <char> DNAbases , vector <int> counts );
vector <char> get_uniques_lc ( vector <char> bases_to_check );

bool contains_chars_not_in_string( std::string test, std::string allowed);

bool is_abba(std::string in);
bool is_baba(std::string in);

void help ();
