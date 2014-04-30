/**************************************************************************************************************************************************************
 *
 * MetaPhase.cc
 *
 * TODO: doc
 *
 *
 * INPUTS:
 * -- The draft assembly (fasta)
 * -- One or more libraries of Hi-C reads (fastq, 2 files per library)
 * -- Optional: The reads that went into the draft assembly (fastq, 2 files per library)
 * -- Optional validation: A set of reference genomes for species known to be in the sample (fasta)
 *
 * POTENTIAL THINGS YOU CAN DO:
 * -- Produce visual output in Cytoscape (MetaAssembly::WriteXGMML)
 *
 *
 * Josh Burton
 * November 2013
 *
 *************************************************************************************************************************************************************/



// C libraries
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// STL declarations
#include <ctime>
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
#include <boost/algorithm/string.hpp> // split
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>


// Modules in ~/include (must add -L~/include and -lJ<module> to link)
#include "TimeMem.h"




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




// Global parameters.  These will eventually be de-hardwired and put in an *.INI file so we don't have to recompile to change them.
// All relative file and directory names are relative to the MetaPhase main directory, ROOT_DIR.
struct RunParams {
  static string SCENARIO; // MY, M2, poplar
  static string ROOT_DIR;
  static string OUTPUT_DIR; // relative to ROOT_DIR
  static string ASSEMBLY_FASTA; // relative to ROOT_DIR
  static string REFS_CSV; // relative to current directory
  static string HIC_LIBS_CSV; // relative to current directory
  static string REFS_DIR; // relative to ROOT_DIR
  static string HIC_DIR; // relative to ROOT_DIR
  static string SHOTGUN_COVERAGE_BAM; // relative to ROOT_DIR
  static string JARVIS_PATRICK_K; // integer
  static string N_CLUSTERS; // integer
  static string MIN_CLUSTER_NORM; // integer
  static string XGMML_FILE; // absolute
  static string LACHESIS_DIR; // relative to OUTPUT_DIR
};

string RunParams::SCENARIO = "M2"; // TEMP: MASTER SWITCH: MY, M2, poplar
string RunParams::ROOT_DIR = "/net/shendure/vol10/jnburton/src/MetaPhase";
string RunParams::OUTPUT_DIR = "out/" + RunParams::SCENARIO;
string RunParams::ASSEMBLY_FASTA = "assembly/MY.ASM2.fasta";
string RunParams::REFS_CSV = "input/" + RunParams::SCENARIO + ".refs.csv";
string RunParams::HIC_LIBS_CSV = "input/" + RunParams::SCENARIO + ".HiC_libs.csv";
string RunParams::REFS_DIR = "refs";
string RunParams::HIC_DIR = "HiC/MY/to_ASM2";
string RunParams::SHOTGUN_COVERAGE_BAM = "assembly/MY/fastq/MetaYeast.to_ASM2.bam";
string RunParams::JARVIS_PATRICK_K = "100";
string RunParams::N_CLUSTERS = "12"; // set to 1 to get the enrichment curve
string RunParams::MIN_CLUSTER_NORM = "25";
string RunParams::XGMML_FILE = "/net/gs/vol2/home/jnburton/public_html/_" + RunParams::SCENARIO + ".xgmml";
string RunParams::LACHESIS_DIR = "Lachesis";





int main( int argc, char * argv[] )
{
  cout << Time() << ": MetaPhase!" << endl;
  
  // Reset options for the other scenarios
  if ( RunParams::SCENARIO == "M2" ) {
    
    bool NEW = true;  // if true, then use the old 18-species assembly, but with new alignments (never use the 19s assembly because R. graminis isn't there!)
    if (NEW) {
      RunParams::OUTPUT_DIR = "out/M2/sim.10k";
      RunParams::ASSEMBLY_FASTA = "assembly/M2/sim.10k.fasta";
      RunParams::HIC_DIR = "HiC/M2/to_sim.10k";
      RunParams::N_CLUSTERS = "18";
      RunParams::MIN_CLUSTER_NORM = "25";
    }
    else {
      RunParams::OUTPUT_DIR = "out/M2/sim_10k";
      RunParams::ASSEMBLY_FASTA = "assembly/M2/sim.10k.fasta";
      RunParams::HIC_DIR = "HiC/M2/old_44bp/to_sim_10k";
      RunParams::N_CLUSTERS = "18";
    }
    
    RunParams::SHOTGUN_COVERAGE_BAM = "";
  }
  else if ( RunParams::SCENARIO == "poplar" ) {
    RunParams::ASSEMBLY_FASTA = "assembly/poplar/PEAR.w_AAGCTT.fasta";
    RunParams::HIC_DIR = "HiC/poplar/to_PEAR.w_AAGCTT";
    RunParams::SHOTGUN_COVERAGE_BAM = "";
    RunParams::N_CLUSTERS = "3"; // TEMP
    RunParams::MIN_CLUSTER_NORM = "1000";
  }
  else assert( RunParams::SCENARIO == "MY" );
  
  // TEMP: other run-time parameters
  const bool load_from_SAM = true; // must be true if do_clustering = true; also necessary for Lachesis output
  const bool do_clustering = false;
  const bool do_bootstrap = false;
  const bool do_merge = false;
  const bool no_output_files = false;
  if ( do_clustering || !no_output_files ) assert( load_from_SAM );
  
  
  // Set up output directory structure.  This includes the main output directory at <ROOT_DIR>/<OUTPUT_DIR>, along with the subdirectory cached_data.
  string out_dir = RunParams::ROOT_DIR + "/" + RunParams::OUTPUT_DIR;
  system( ( "mkdir -p " + out_dir ).c_str() );
  string cache_dir = out_dir + "/cached_data";
  system( ( "mkdir -p " + cache_dir ).c_str() );
  
  
  
  
  string assembly_fasta = RunParams::ROOT_DIR + "/" + RunParams::ASSEMBLY_FASTA;
  string shotgun_BAM = RunParams::ROOT_DIR + "/" + RunParams::SHOTGUN_COVERAGE_BAM;
  if ( RunParams::SHOTGUN_COVERAGE_BAM == "" ) shotgun_BAM = "";
  
  
  
  // Create the TrueMapping object, and load information about the references into it.
  TrueMapping truth( assembly_fasta, RunParams::REFS_CSV, RunParams::ROOT_DIR + "/" + RunParams::REFS_DIR, cache_dir + "/TrueMapping.txt" );
  
  
  // Print TargetCoverages.  This is useful in creating big ol' histograms of read coverages for supp. figures (see make_histos.sh)
  if (0) {
    truth.ReportTargetCoverages( RunParams::ROOT_DIR + "/assembly/MY/fastq/MetaYeast.to_ASM2.bam" );
    return 0;
  }
  
  
  
  
  
  // Get the library names from the libraries CSV file.
  vector<struct HiCLib> HiC_libs = ParseHiCLibsCSV( RunParams::HIC_LIBS_CSV, RunParams::ROOT_DIR + "/" + RunParams::HIC_DIR );
  
  
  int N_clusters = boost::lexical_cast<int>( RunParams::N_CLUSTERS );
  
  // Run the clustering algorithm on each library in turn.
  vector<ClusteringResult> clustering_results;
  for ( size_t i = 0; i < HiC_libs.size(); i++ ) {
    string lib_name = HiC_libs[i].name;
    string RE_site  = HiC_libs[i].RE_site;
    string SAM_file = HiC_libs[i].SAM_file;
    
    
    // Create a MetaAssembly object from this assembly.
    MetaAssembly ma( assembly_fasta, RE_site );
    ma.LoadTrueMapping( &truth );
    if ( load_from_SAM )
      ma.LoadFromSAM( SAM_file );
    
    
    // Truth-based pre-analysis.
    //ma.FindClusterDensities();
    

    // Run clustering (or just read the clustering result from file.)
    string clusters_file = cache_dir + "/clusters." + lib_name + ".K" + RunParams::N_CLUSTERS + ".txt";
    
    
    if ( do_clustering ) {
      if ( do_bootstrap ) ma.Bootstrap();
      ma.PrepLinkMatrix( boost::lexical_cast<int>(RunParams::JARVIS_PATRICK_K), shotgun_BAM );
      ma.Cluster( N_clusters, boost::lexical_cast<int>(RunParams::MIN_CLUSTER_NORM) );
      if ( !no_output_files ) {
	cout << Time() << ": Writing to clusters file " << clusters_file << endl;
	ma._clusters.WriteFile( clusters_file );
      }
    }
    else {
      cout << Time() << ": Reading from clusters file " << clusters_file << endl;
      ma._clusters.ReadFile( clusters_file );
    }
    
    clustering_results.push_back( ma._clusters );
    
    ma.FindClusterDensities(false);
    
    
    if ( load_from_SAM ) {
      //ma.FindIntraClusterEnrichment();
      //ma.FindMultiSpeciesContigs(); // find and also evaluate multi-species contigs
    }
    
    // Report on the clustering result.
    ma._clusters.ReorderClustersByRefs( truth );
    ma._clusters.Print( truth );
    if ( load_from_SAM ) ma.ReportUnclusteredContigs();
    
    if ( no_output_files ) continue;
    
    // Output files for certain clusters, to be used by Lachesis.
    if ( RunParams::SCENARIO == "poplar" ) { // poplar: all clusters, for tetramer testing
      for ( int i = 0; i < ma._clusters.NClusters(); i++ )
	ma.WriteClusterFasta( out_dir + "/clusters." + boost::lexical_cast<string>(i) + ".fasta", i, false );
    }
    
    else if ( RunParams::SCENARIO == "M2" ) { // Menagerie2: M. maripaludis, B. subtilis
      ma.WriteClusterFasta( out_dir + "/" + RunParams::LACHESIS_DIR + ".MM/contigs.fasta", 8,  true ); // 8 = M. maripaludis
      ma.WriteClusterFasta( out_dir + "/" + RunParams::LACHESIS_DIR + ".VF/contigs.fasta", 10, true ); // 10 = V. fischeri
      ma.WriteClusterFasta( out_dir + "/" + RunParams::LACHESIS_DIR + ".RP/contigs.fasta", 15, true ); // 15 = R. palustris
      ma.WriteClusterFasta( out_dir + "/" + RunParams::LACHESIS_DIR + ".BS/contigs.fasta", 17, true ); // 17 = B. subtilis
    }
    else { // MetaYeast case
      ma.WriteClusterFasta( out_dir + "/" + RunParams::LACHESIS_DIR + ".KW/contigs.fasta", 12, true ); // 12 = K. wickerhamii
      ma.WriteClusterFasta( out_dir + "/" + RunParams::LACHESIS_DIR + ".SS/contigs.fasta", 14, true ); // 14 = S. stipitis
      
      ma.WriteClusterFasta( out_dir + "/unclustered.fasta", -1, false ); // write unclustered contigs
    }
    
    // Make plots.
    if ( !load_from_SAM ) continue;
    
    // Draw Figure 2a (the big clustering image.)
    string Fig2a_file = "~/public_html/graph." + RunParams::SCENARIO + ".png";
    ma.DrawFigure2a( Fig2a_file );
    
    ma._clusters.DrawFigure2bc( truth, RunParams::SCENARIO == "MY" );
    if ( RunParams::SCENARIO == "MY" ) {
      system( ( "mv ~/public_html/Fig2b.jpg ~/public_html/Fig2b." + lib_name + ".jpg" ).c_str() );
      system( ( "mv ~/public_html/Fig2c.jpg ~/public_html/Fig2c." + lib_name + ".jpg" ).c_str() );
    }
    ma.WriteXGMML( RunParams::XGMML_FILE, 2 );
    
    
    
  }
  
  
  if ( !do_merge ) { cout << Time() << ": Done!" << endl; return 0; }
  
  // Merge the two results.
  assert( clustering_results.size() == 2 );
  ClusteringResult merged_clusters = MergeClusteringResults1( clustering_results[0], clustering_results[1] );
  
  // Merge the two results (TEMP: new method)
  merged_clusters = MergeClusteringResults2( clustering_results, N_clusters );
  
  // Report on the merged result.
  merged_clusters.ReorderClustersByRefs( truth );
  merged_clusters.Print( truth );
  
  if ( no_output_files ) { cout << Time() << ": Done!" << endl; return 0; }
  
  //merged_clusters.DrawChart( truth, false, true, RunParams::SCENARIO == "MY", 500000 );
  //merged_clusters.DrawChart( truth, true,  true, RunParams::SCENARIO == "MY", 500000 );
  //merged_clusters.DrawFigure2a();
  merged_clusters.DrawFigure2bc( truth, RunParams::SCENARIO == "MY" );
  system( "mv ~/public_html/Fig2b.jpg ~/public_html/Fig2b.merged.jpg" );
  system( "mv ~/public_html/Fig2c.jpg ~/public_html/Fig2c.merged.jpg" );
  
  // Write the merged result to file.
  string clusters_file = cache_dir + "/clusters.HN.K" + boost::lexical_cast<string>( merged_clusters.NClusters() ) + ".txt";
  merged_clusters.WriteFile( clusters_file );
  
  
  // Create a MetaAssembly from this merged result, just to write output files to be input into Lachesis.
  MetaAssembly ma( assembly_fasta, "AAGCTT_CCATGG" );
  ma.SetClusters( merged_clusters );
  ma.LoadTrueMapping( &truth );
  
  ma.WriteClusterFasta( out_dir + "/" + RunParams::LACHESIS_DIR + "/contigs.fasta", 14, true ); // 14 = S. stipitis
  
  //ma.LoadFromSAM( RunParams::ROOT_DIR + "/" + RunParams::HIC_DIR +  "/Nco1.all_lens.REduced.500.bam" );
  //ma.LoadFromSAM( RunParams::ROOT_DIR + "/" + RunParams::HIC_DIR + "/Hind3.all_lens.REduced.500.bam" );
  //ma.FindMultiSpeciesContigs(); // find and also evaluate multi-species contigs (needs link data via LoadFromSAM())
  
  
  cout << Time() << ": Done!" << endl;
  return 0;
}



