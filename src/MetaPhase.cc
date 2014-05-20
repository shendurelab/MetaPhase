/**************************************************************************************************************************************************************
 *
 * MetaPhase.cc
 *
 * This is the top-level module of the MetaPhase software program.  This is where files are loaded in and objects are created.
 * The other modules (each of which have their own .cc and .h file) are as follows:
 * MetaAssembly:           The main object that contains Hi-C link data (in the form of a LinkMatrix), on which clustering and other algorithms are performed.
 * LinkMatrix:             A glorified Boost::ublas matrix object; an upper triangular matrix containing (possibly normalized) Hi-C link density among contigs.
 * HierarchicalClustering: The agglomerative hierarchical clustering algorithm.
 * ClusteringResult:       A putative assignment of contigs to clusters, produced by the clustering algorithm.
 * TrueMapping:            A record of which contigs in the draft assembly map to which reference genomes.
 * TextFileParsers:        Several subroutines for parsing input text files and converting the results to native C++ data types.
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
#include <iostream>
#include <fstream>
using namespace std;

// Local includes
#include "MetaAssembly.h"
#include "HiCLib.h"
#include "TrueMapping.h"

// Boost includes
#include <boost/lexical_cast.hpp>


// Modules in ~/include (must add -L~/include and -lJ<module> to link)
#include "TimeMem.h"







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

string RunParams::SCENARIO = "MY"; // TEMP: MASTER SWITCH: MY, M2, poplar
string RunParams::ROOT_DIR = "/net/shendure/vol10/jnburton/src/MetaPhase";
string RunParams::OUTPUT_DIR = "out/" + RunParams::SCENARIO;
string RunParams::ASSEMBLY_FASTA = "assembly/MY.ASM2.fasta";
string RunParams::REFS_CSV = "input/" + RunParams::SCENARIO + ".refs.csv";
string RunParams::HIC_LIBS_CSV = "input/" + RunParams::SCENARIO + ".HiC_libs.csv";
string RunParams::REFS_DIR = "refs";
string RunParams::HIC_DIR = "HiC/MY/to_ASM2";
string RunParams::SHOTGUN_COVERAGE_BAM = "assembly/MY/fastq/MetaYeast.to_ASM2.bam";
string RunParams::JARVIS_PATRICK_K = "100";
string RunParams::N_CLUSTERS = "12"; // set to 1 to get the enrichment curve! - to get enrichments from stdout: grep "E\(" c | awk '{print $8,$5}' | uniq -f1 | awk '{print $2,$1}' | tail -n50
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
  const bool load_from_SAM = true; // necessary for the clustering algorithm and many output files
  const bool do_clustering = false; // if this is true, load_from_SAM must be true
  const bool do_bootstrap = false;
  const bool do_merge = false;
  // TEMP: output files - set to all false to make no output files
  const bool output_clusters = false; // make file w/ contig clustering assignments; has no effect unless do_clustering = true
  const bool output_network_image = true; // make the big igraph image of the clustering result - time-consuming!
  const bool output_heatmaps = true; // make heatmaps indicating true contig placements
  const bool output_cluster_fastas = false; // fastas of individual clusters - used by Lachesis
  
  // Enforce sanity in run-time parameters.
  if ( output_network_image || output_heatmaps || output_cluster_fastas ) assert( load_from_SAM );
  
  
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
      if ( output_clusters ) {
	cout << Time() << ": Writing to clusters file " << clusters_file << endl;
	ma._clusters.WriteFile( clusters_file );
      }
    }
    else {
      cout << Time() << ": Reading from clusters file " << clusters_file << endl;
      ma._clusters.ReadFile( clusters_file );
    }
    
    clustering_results.push_back( ma._clusters );
    
    // Optional (runtime-consuming) report functions.
    if ( 0 ) {
      //ma.FindClusterDensities(false);
      if ( load_from_SAM ) {
	//ma.FindIntraClusterEnrichment();
	//ma.FindMultiSpeciesContigs(); // find and also evaluate multi-species contigs
	//ma.ReportUnclusteredContigs();
      }
    }
    
    
    // Report on the clustering result.
    ma._clusters.ReorderClustersByRefs( truth );
    ma._clusters.Print( truth );
    
    // Make plots.
    if ( !load_from_SAM ) continue;
    
    // Draw the big image of the clustering network.
    if ( output_network_image ) ma.DrawClusterNetwork( "~/public_html/graph." + RunParams::SCENARIO + ".png", RunParams::SCENARIO );
    
    // Draw the heatmaps indicating which clusters contain contigs truly belonging to which species.
    if ( output_heatmaps ) {
      ma._clusters.DrawTruthHeatmaps( truth, RunParams::SCENARIO == "MY" );
      if ( RunParams::SCENARIO == "MY" ) {
	system( ( "mv ~/public_html/Fig2b.jpg ~/public_html/Fig2b." + lib_name + ".jpg" ).c_str() );
	system( ( "mv ~/public_html/Fig2c.jpg ~/public_html/Fig2c." + lib_name + ".jpg" ).c_str() );
      }
    }
    
    // Output files for certain clusters, to be used by Lachesis.
    if ( output_cluster_fastas ) {
      
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
      
    }
    
    // XGMML file for Cytoscape: I don't use this anymore
    //ma.WriteXGMML( RunParams::XGMML_FILE, 2 );
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
  
  if ( output_heatmaps ) {
    merged_clusters.DrawTruthHeatmaps( truth, RunParams::SCENARIO == "MY" );
    system( "mv ~/public_html/Fig2b.jpg ~/public_html/Fig2b.merged.jpg" );
    system( "mv ~/public_html/Fig2c.jpg ~/public_html/Fig2c.merged.jpg" );
  }
  
  // Write the merged result to file.
  if ( output_clusters ) {
    string clusters_file = cache_dir + "/clusters.HN.K" + boost::lexical_cast<string>( merged_clusters.NClusters() ) + ".txt";
    merged_clusters.WriteFile( clusters_file );
  }
  
  
  // Create a MetaAssembly from this merged result, just for output purposes.
  MetaAssembly ma( assembly_fasta, "AAGCTT_CCATGG" );
  ma.SetClusters( merged_clusters );
  ma.LoadTrueMapping( &truth );
  ma.LoadFromSAM( RunParams::ROOT_DIR + "/" + RunParams::HIC_DIR +  "/Nco1.all_lens.REduced.500.bam" );
  ma.LoadFromSAM( RunParams::ROOT_DIR + "/" + RunParams::HIC_DIR + "/Hind3.all_lens.REduced.500.bam" );
  ma.FindMultiSpeciesContigs(); // find and also evaluate multi-species contigs (needs link data via LoadFromSAM())
  if ( output_network_image ) ma.DrawClusterNetwork( "~/public_html/graph." + RunParams::SCENARIO + ".merged.png", RunParams::SCENARIO );
  
  
  cout << Time() << ": Done!" << endl;
  return 0;
}



