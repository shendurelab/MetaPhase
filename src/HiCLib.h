/**************************************************************************************************************************************************************
 *
 * HiCLib.h
 *
 * HiCLib: A simple struct containing info about a Hi-C library - all the same info as in one line of the HiC_libs CSV.
 *
 *
 *************************************************************************************************************************************************************/


#ifndef _HIC_LIB__H
#define _HIC_LIB__H


// C libraries
#include <stdio.h>
#include <assert.h>

// STL declarations
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
using namespace std;

// Local includes
#include "MetaAssembly.h"
#include "TrueMapping.h"
#include "TextFileParsers.h" // TokenizeCSV

// Boost includes
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>





// HiCLib: A simple struct containing info about a Hi-C library - all the same info as in one line of the HiC_libs CSV.
// This struct is filled in the function ParseHiCLibsCSV.
struct HiCLib {
  string name; // used only in the filename of the clustering cache file
  string RE_site; // the restriction enzyme site sequence
  string SAM_file; // absolute location of a SAM/BAM file containing this library aligned to the draft assembly
};



// ParseHiCLibsCSV: Parse a HiC_libs CSV file and produce a vector<HiCLib>.
// This function includes tons of sanity checks and verbose output, so that any kind of user error in creating the CSV file will be immediately reported.
vector<HiCLib>
ParseHiCLibsCSV( const string & CSV_file, const string & HiC_dir )
{
  vector<HiCLib> libs;
  set<string> names;
  
  if ( !boost::filesystem::is_regular_file( CSV_file ) ) {
    cout << "ERROR: Can't find input CSV file `" << CSV_file << "'" << endl;
    exit(1);
  }
  
  vector< vector<string> > csv_tokens;
  TokenizeCSV( CSV_file, csv_tokens );
  for ( size_t i = 0; i < csv_tokens.size(); i++ ) {
    vector<string> & tokens = csv_tokens[i];
    
    // Skip commented lines.
    assert( tokens.size() > 0 );
    assert( tokens[0].size() > 0 );
    if ( tokens[0][0] == '#' ) continue;
    
    
    // Do a bunch of sanity checks.
    string error_str = "ERROR: Input CSV file `" + CSV_file + "', line " + boost::lexical_cast<string>(i);
    
    if ( tokens.size() != 3 ) {
      cout << error_str << ": Should see 3 tokens, instead saw " << tokens.size() << "." << endl;
      exit(1);
    }
    
    if ( find( names.begin(), names.end(), tokens[0] ) != names.end() ) {
      cout << error_str << ": Abbreviated name `" << tokens[0] << "' is already used by another line in this file.  These names should be unique!" << endl;
      exit(1);
    }
    
    if ( tokens[1].find_first_not_of( "ACGTN" ) != string::npos ) {
      cout << error_str << ": Restriction enzyme site `" << tokens[1] << "' should not contain any characters other than A,C,G,T,N!" << endl;
      exit(1);
    }
    
    string SAM = HiC_dir + "/" + tokens[2];
    if ( !boost::filesystem::is_regular_file( SAM ) ) {
      cout << error_str << ": Can't find SAM file `" << SAM << "'." << endl;
      exit(1);
    }
    
    // Make a HiCLib object out of these strings.
    HiCLib lib;
    lib.name = tokens[0];
    lib.RE_site = tokens[1];
    lib.SAM_file = SAM;
    libs.push_back(lib);
    
    // Keep track of library names that have been seen, so that the same name won't be used twice.
    names.insert( tokens[0] );
  }
  
  
  return libs;
}




#endif
