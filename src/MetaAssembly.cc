// For documentation, see MetaAssembly.h
#include "MetaAssembly.h"
#include "TrueMapping.h"
#include "TextFileParsers.h" // TokenizeFile, GetFastaNames
#include "LinkMatrix.h" // LinkMatrixInt
#include "HierarchicalClustering.h" // AgglomerativeHierarchicalClustering
#include "DensityClustering.h" // DensityClustering

#include <assert.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm> // fill_n, count

// Boost includes
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

// Modules in ~/include (must add -L~/include and -lJ<module> to link)
#include "TimeMem.h"
#include "gtools/N50.h"
#include "gtools/SAMStepper.h" // SAMStepper, NTargetsInSAM






MetaAssembly::MetaAssembly( const string & assembly_fasta, const string & RE_site_seq )
  : _assembly_fasta( assembly_fasta ),
    _truth( NULL ),
    _RE_sites_file( assembly_fasta + ".counts_" + RE_site_seq + ".txt" )
{
  cout << Time() << ": Creating a MetaAssembly from " << assembly_fasta << endl;
  assert( boost::filesystem::is_regular_file( assembly_fasta ) );
  
  // Get the names of the assembly contigs from the fasta file.
  _contig_names = GetFastaNames( assembly_fasta );
  _N_contigs = _contig_names.size();
  PRINT( _N_contigs );
  assert( _N_contigs > 0 );
  
  
  // Create the FastaSize file, if necessary.
  string FastaSize_file = assembly_fasta + ".FastaSize";
  if ( !boost::filesystem::is_regular_file( FastaSize_file ) ) {
    string cmd = "FastaSize " + assembly_fasta + " > " + FastaSize_file;
    cout << Time() << ": " << cmd << endl;
    system( cmd.c_str() );
    assert( boost::filesystem::is_regular_file( FastaSize_file ) );
  }
  
  _contig_lengths.resize( _N_contigs );
  
  vector< vector<string> > FastaSize_tokens;
  TokenizeFile( assembly_fasta + ".FastaSize", FastaSize_tokens, true );
  assert( _N_contigs == (int) FastaSize_tokens.size() - 1 ); // the FastaSize file has one line for each contig, plus a summary line at the end
  
  // Get the contig lengths from the fasta.FastaSize file.  For each line in the FastaSize file (before the last), there should be 3 tokens: a blank one
  // (whitespace), a number indicating length, and a chromosome name.
  for ( int i = 0; i < _N_contigs; i++ ) {
    vector<string> & t = FastaSize_tokens[i];
    assert( t.size() == 3 );
    assert( t[0] == "" );
    assert( t[2] == _contig_names[i] );
    _contig_lengths[i] = boost::lexical_cast<int>( t[1] );
  }
  
  
  
  // Get the number of RE sites per contig.  A file may or may not have already been generated with this information (via CountMotifsInFasta.pl.)
  // If it hasn't, we need to run CountMotifsInFasta now.
  if ( !boost::filesystem::is_regular_file( _RE_sites_file ) ) {
    string cmd = "CountMotifsInFasta.pl " + assembly_fasta + " " + RE_site_seq;
    cout << Time() << ": To create file " << _RE_sites_file << ": " << cmd << endl;
    system( cmd.c_str() );
    assert( boost::filesystem::is_regular_file( _RE_sites_file ) ); // if this fails, the CountMotifsInFasta.pl script didn't run correctly
  }
  
  
  // Get token 1 (zero-indexed) from each line of the RE sites file.
  _contig_RE_sites = ParseTabDelimFile<int>( _RE_sites_file, 1 );
  assert( (int) _contig_RE_sites.size() == _N_contigs );
  
  // Modify contig_RE_sites so there are no zeros to divide by.
  for ( int i = 0; i < _N_contigs; i++ )
    _contig_RE_sites[i]++;
  
  
  _SAM_files.clear();
}









// Load a SAM/BAM file and put its links into this MetaAssembly.
void
MetaAssembly::LoadFromSAM( const string & SAM_file )
{
  cout << Time() << ": Loading from SAM file " << SAM_file << endl;
  _SAM_files.push_back( SAM_file );
  
  if ( NTargetsInSAM( SAM_file ) != _N_contigs ) {
    cout << "ERROR: SAM file " << SAM_file << " has " << NTargetsInSAM( SAM_file ) << " target contigs, but genome assembly fasta has " << _N_contigs << " contigs.  Are your files mismatched?" << endl;
    exit(1);
  }
  
  bool verbose = true;
  
  // If this is the first SAM file, set up the link matrix to handle the data as it comes in.  The matrix is upper triangular, so matrix(i,j) = 0 unless j > i.
  // It's a sparse matrix so it doesn't take up too much memory.
  if ( _link_matrix.size1() == 0 ) _link_matrix.Resize( _N_contigs );
  
  int N_pairs_here = 0, N_links_here = 0;
  vector<bool> seen_link( _N_contigs, false );
  
  // Set up a SAMStepper object to read in the alignments.
  SAMStepper stepper(SAM_file);
  stepper.FilterAlignedPairs(); // Only look at read pairs where both reads aligned to the draft assembly.
  
  
  // Loop over all pairs of alignments in the SAM file.
  // Note that the next_pair() function assumes that all reads in a SAM file are paired, and the two reads in a pair occur in consecutive order.
  for ( pair< bam1_t *, bam1_t *> aligns = stepper.next_pair(); aligns.first != NULL; aligns = stepper.next_pair() ) {
    
    if ( verbose && stepper.N_aligns_read() % 1000000 == 0 ) cout << "." << flush;
    
    const bam1_core_t & c1 = aligns.first->core;
    const bam1_core_t & c2 = aligns.second->core;
    
    // If the two reads align to the same contig, the link isn't informative, so skip it.
    if ( c1.tid == c2.tid ) continue;
    
    // Ignore reads with mapping quality 0.
    if ( c1.qual == 0 || c2.qual == 0 ) continue;
    
    // Sanity checks to make sure the read pairs appear as they should in the SAM file.  If they don't, the read pair is corrupt, so skip it.
    if ( c1.tid != c2.mtid ) continue;
    if ( c2.tid != c1.mtid ) continue;
    if ( c1.pos != c2.mpos ) continue;
    if ( c2.pos != c1.mpos ) continue;
    if ( c1.tid >= _N_contigs ) continue;
    if ( c2.tid >= _N_contigs ) continue;
    
    // Rearrange the order of the contig IDs if necessary, so that ID1 < ID2.
    int ID1, ID2;
    if ( c1.tid < c2.tid ) { ID1 = c1.tid; ID2 = c2.tid; }
    else                   { ID1 = c2.tid; ID2 = c1.tid; }
    //PRINT2( ID1, ID2 );
    
    
    // Add to tallies.
    if ( _link_matrix(ID1,ID2) == 0 ) N_links_here++;
    N_pairs_here++;
    seen_link[ID1] = true;
    seen_link[ID2] = true;
    
    // Add the data to the matrix.
    _link_matrix(ID1,ID2) += 1;
  }
  
  
  if ( verbose ) cout << endl;
  
  _link_matrix.N_links += N_links_here;
  _link_matrix.N_pairs += N_pairs_here;
  
  cout << Time() << ": N aligns/pairs read from this file: " << stepper.N_aligns_read() << "/" << stepper.N_pairs_read() << "; N pairs used: " << N_pairs_here << "; N contigs = " << _N_contigs << "; matrix now has " << _link_matrix.N_links << " elements with sum = " << _link_matrix.N_pairs << endl;
  
  
  int N_singletons = count( seen_link.begin(), seen_link.end(), false );
  PRINT( N_singletons );
}





// Find the average link density within clusters.
// If use_truth = true, this a truth-based pre-analysis: it finds for each species the average linkage between contigs that align uniquely to that species.
// If use_truth = false, this is a non-truth-based post-analysis, operating on the derived clusters.
void
MetaAssembly::FindClusterDensities( const bool use_truth ) const
{
  cout << Time() << ": FindClusterDensities!" << endl;
  
  assert( _link_matrix.has_data() ); // if this fails, you need to load data into the MetaAssembly via LoadFromSAM()
  if ( use_truth ) assert( _truth != NULL ); // if this fails, call LoadTrueMapping()
  else assert( !_clusters.empty() ); // if this fails, call Cluster()
  
  int N_clusters = use_truth ? _truth->NRefs() : _clusters.NClusters();
  PRINT( N_clusters );
  
  vector<int> contig_to_cluster( _N_contigs, -1 );
  vector<int> cluster_size( N_clusters, 0 );
  vector<int> cluster_len( N_clusters, 0 ); // "length" in RE sites, not bp
  
  
  if ( use_truth )
    
    // Loop over all reference genome IDs and find contigs mapping uniquely to that reference.
    for ( int ref_ID = 0; ref_ID < N_clusters; ref_ID++ ) {
      
      // Find all contigs that align uniquely to this reference.
      vector<int> unique_cIDs;
      for ( int i = 0; i < _N_contigs; i++ )
	if ( _truth->QOnRefOnly( i, ref_ID ) )
	  contig_to_cluster.at(i) = ref_ID;
    }
  else
    
    for ( int i = 0; i < _N_contigs; i++ )
      contig_to_cluster[i] = _clusters.cluster_ID(i);
  
  
  // Fill the cluster_size and cluster_len vectors.
  for ( int i = 0; i < _N_contigs; i++ ) {
    int cID = contig_to_cluster[i];
    if ( cID == -1 ) continue;
    cluster_size[cID]++;
    cluster_len [cID] += _contig_RE_sites[i];
  }
  
  
  // Now loop over the matrix and find all intra-cluster links.
  vector<int> N_links( N_clusters, 0 );
  
  for ( LinkMatrixInt::const_iterator1 it1 = _link_matrix.begin1(); it1 != _link_matrix.end1(); ++it1 )
    for ( LinkMatrixInt::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2 ) {
      int cluster1 = contig_to_cluster[ it2.index1() ];
      int cluster2 = contig_to_cluster[ it2.index2() ];
      if ( cluster1 != -1 && cluster1 == cluster2 )
	N_links[cluster1] += *it2;
    }
  
  
  // Normalize the link numbers by the contig lengths.
  vector<double> N_links_avg;
  for ( int cID = 0; cID < N_clusters; cID++ )
    N_links_avg.push_back( double( N_links[cID] ) / cluster_len[cID] / cluster_len[cID] );
  
  // Rescale the link density numbers by the magnitude of the average number, to make the numbers more human-readable.
  vector<double> nla_sort = N_links_avg;
  sort( nla_sort.begin(), nla_sort.end() );
  int log_norm = log10( nla_sort[N_clusters/2] ) - 1;
  PRINT( log_norm );
  for ( int cID = 0; cID < N_clusters; cID++ )
    N_links_avg[cID] /= pow( 10, log_norm );
  
  
  // Print!
  cout << Time() << ": Report!" << endl;
  if ( use_truth )
    for ( int ref_ID = 0; ref_ID < N_clusters; ref_ID++ )
      cout << "Ref #" << ref_ID << " (" << _truth->RefName(ref_ID) << ")\tAvg. linkage = " << N_links_avg[ref_ID] << "\t\tN unique contigs = " << cluster_size[ref_ID] << ", len = " << cluster_len[ref_ID] << endl;
  else
    for ( int cID = 0; cID < N_clusters; cID++ )
      cout << "Cluster #" << cID << "\tPlurality ref: " << _clusters.PluralityRefName( *_truth, cID ) << "\tAvg. linkage = " << N_links_avg[cID] << "\t\tN contigs = " << cluster_size[cID] << ", len = " << cluster_len[cID] << endl;
  
} 




// Prepare the link matrix for clustering, by normalization etc.  This uses _link_matrix and fills _normed_link_matrix.
void
MetaAssembly::PrepLinkMatrix( const int K, const string & shotgun_BAM )
{
  assert( _link_matrix.has_data() ); // if this fails, you need to load data into the MetaAssembly via LoadFromSAM()
  
  static const int ISOLATED_COMPONENT_SIZE = 100; // contigs in components of this size or smaller are ignored (DO NOT SET to 0)
  
  _link_matrix.IsolateNormalizeLimit( _contig_RE_sites, ISOLATED_COMPONENT_SIZE, K, shotgun_BAM, _normed_link_matrix );
  assert( _normed_link_matrix._isolated_contigs.size() > 0 );
}



void
MetaAssembly::Cluster( const int N_clusters, const int MIN_CLUSTER_NORM )
{
  assert( _normed_link_matrix.has_data() ); // if this fails, you need to load data into the MetaAssembly via LoadFromSAM() and then call PrepLinkMatrix()
  assert( _normed_link_matrix._isolated_contigs.size() > 0 );
  
  vector<int> contig_clusters = AgglomerativeHierarchicalClustering( _link_matrix, _normed_link_matrix, _contig_RE_sites, N_clusters, MIN_CLUSTER_NORM );
  
  _clusters = ClusteringResult( contig_clusters, _contig_lengths );
}








/* FindMultiSpeciesContigs: Find contigs that could plausibly belong to more than one cluster, implying that their sequences appear in multiple species.
 *
 * This is a crucial part of the MetaPhase method.  Different species, and especially different strains of the same species, are going to have long stretches
 * of orthologous sequence.  Metagenome assemblers will tend to collapse these into single contigs because they don't have super-long-range contiguity
 * information and cannot make any read-depth-based assumptions about copy number.  Let's call these "multi-species contigs".
 * 
 * Multi-species contigs pose a problem for metagenomic deconvolution.  They will have biologically genuine Hi-C links to other contigs from all of the species
 * to which they belong, and will thus appear to belong to multiple clusters in the MetaAssembly link graph.  Fortunately, there aren't too many multi-species
 * contigs, and their ambiguous signal should be swamped by the unambiguous signal of single-species contigs, so that the Cluster() function can still produce
 * correct clusters.  After clustering is done, non-unique contigs can be identified ex post facto because they should have a linkage signal with one ot more
 * other clusters that is as strong as their intra-cluster signal.  That's the criterion we use here.
 *
 * With the current (sparse) link density provided by the preliminary dataset, our criterion is as follows.  Let C be a contig in cluster A.  We call C a
 * multi-species contig, in cluster B as well as A, if the link density between C and the contigs in cluster B is at least 1 links per megabase of sequence
 * in cluster B.  We would like to normalize by the length of contig C as well, but unfortunately there is not enough signal in the sample, so doing this
 * normalization makes the resulting calls far worse (this has been tested for the SC/SP multi-species contig calls.)
 *
 * This function also includes truth-based evaluation of the multi-species contig calls.
 *
 */
void
MetaAssembly::FindMultiSpeciesContigs()
{
  cout << Time() << ": FindMultiSpeciesContigs" << endl;
  
  assert( !_clusters.empty() );
  assert( _truth != NULL ); // if this fails, call LoadTrueMapping()
  
  const bool USE_RES = true; // if true, SCALING_FACTOR is counted in RE sites, not bp, and all normalization is done by RE sites
  const double SCALING_FACTOR = 1.0e2; // i.e., 100 RE sites (if USE_RES = true) or 100 bp (otherwise); this makes the link densities more human-readable
  const double DENSITY_THRESHOLD = 1; // a contig can be called as multi-species if its link density to a cluster is this times the intra-cluster density
  
  int N_clusters = _clusters.NClusters();
  
  // Keep track of the number of true/false positives, which are determined by comparing the multi-species contig calls to true reference alignments.
  vector< vector<int> > N_true_positives ( N_clusters, vector<int>( N_clusters, 0 ) );
  vector< vector<int> > N_false_positives( N_clusters, vector<int>( N_clusters, 0 ) );
  int N_true_positives_total = 0;
  int N_false_positives_total = 0;
  int N_error_fixes_total = 0;
  int total_len = 0;
  int true_len = 0;
  
  
  // For the SC-SP case specifically (refs 0, 4) keep track of which contigs are called as multi-species.  This will let us find the false negative rate later.
  int N_refs_04_called = 0;
  int len_refs_04_called = 0;
  
  // As a pre-calculation step, find the size and the Hi-C link density (per SCALING_FACTOR) within each cluster.
  vector<int>    cluster_N_links( N_clusters, 0 ); // number of links within each cluster
  vector<int>    cluster_norms  ( N_clusters, 0 ); // length of each cluster (in either bp or RE sites, depending on USE_RES)
  vector<int>    cluster_norm2s ( N_clusters, 0 );
  vector<double> cluster_density( N_clusters, 0 ); // this is the only output we need
  for ( int i = 0; i < _N_contigs; i++ ) {
    int cID = _clusters.cluster_ID(i);
    if ( cID == -1 ) continue; // contig wasn't clustered
    
    // Add this contig to its cluster norms.
    cluster_norms[cID] += ( USE_RES ? _contig_RE_sites[i] : _contig_lengths[i] );
    
    // Find all other contigs in this cluster, and tally up the links between them.
    for ( int j = i+1; j < _N_contigs; j++ )
      if ( _clusters.cluster_ID(j) == cID )
	cluster_N_links[cID] += _link_matrix(i,j);
  }
  
  // Calculate the total range ("norm^2") over which clusters can have inter-contig links.  This is non-trivial, because intra-contig links are not counted.
  // The range is N^2 / sum( n^2 ), where N = the total cluster norm and n = the norms of the individual contigs in the cluter.
  for ( int i = 0; i < N_clusters; i++ ) {
    cluster_norm2s[i] = cluster_norms[i] * cluster_norms[i];
    assert( cluster_norm2s[i] > 0 ); // check for integer overflow
  }
  
  for ( int i = 0; i < _N_contigs; i++ ) {
    int cID = _clusters.cluster_ID(i);
    if ( cID == -1 ) continue; // contig wasn't clustered
    int norm = ( USE_RES ? _contig_RE_sites[i] : _contig_lengths[i] );
    cluster_norm2s[cID] -= norm * norm;
  }
  
  
  // Now calculate the Hi-C link densities.
  for ( int i = 0; i < N_clusters; i++ ) {
    cluster_density[i] = cluster_N_links[i] * SCALING_FACTOR * SCALING_FACTOR / cluster_norm2s[i];
    PRINT4( _clusters.ClusterSize(i), cluster_norm2s[i], cluster_N_links[i], cluster_density[i] );
  }
  
  
  
  
  
  // Main loop: Loop over all contigs.
  for ( int i = 0; i < _N_contigs; i++ ) {
    int cID = _clusters.cluster_ID(i);
    if ( cID == -1 ) continue; // contig wasn't clustered; don't attempt to cluster it here, that's not what this function is for
    
    
    
    bool seen_inter_link = false;
    
    // For this contig, determine how many links it has to contigs in every cluster, including its own.
    vector<int> N_links_to_cluster( N_clusters, 0 );
    
    // Loop over all of the non-zero entries in the link_matrix row/column corresponding to this contig.
    // TODO: this may be slow, but I can't figure out a better (iterator-based) method that doesn't have a huge memory footprint
    for ( int j = 0; j < _N_contigs; j++ ) {
      int cID2 = _clusters.cluster_ID(j);
      if ( cID2 == -1 ) continue;
      int N_links = _link_matrix( min(i,j), max(i,j) );
      if ( N_links == 0 ) continue;
      
      N_links_to_cluster[cID2] += N_links;
      if ( cID != cID2 ) seen_inter_link = true;
    }
    
    
    if ( !seen_inter_link ) continue; // this contig doesn't have any links outside its cluster, so we have no evidence that it's a multi-species contig
    
    
    // For each cluster, consider the possibility that this contig may belong in this cluster (in addition to its own.)
    vector<double> link_density_to_cluster( N_clusters, 0 );
    //PRINT3( i, _contig_lengths[i], cID );
    for ( int j = 0; j < N_clusters; j++ ) {
      if ( j == cID ) continue; // no need to look for a contig's connection to its own cluster!
      if ( _clusters.ClusterSize(j) < 10 ) continue; // no need to look for connections to tiny clusters, which don't represent the plurality of a species
      
      // Find the normalized density of Hi-C links between this contig and this cluster.
      double link_density_to_cluster = double( N_links_to_cluster[j] );
      link_density_to_cluster *= SCALING_FACTOR / cluster_norms[j];
      link_density_to_cluster *= SCALING_FACTOR / ( USE_RES ? _contig_RE_sites[i] : _contig_lengths[i] );
      //PRINT3( link_density_to_cluster, cluster_density[j], _clusters.ClusterSize(j) );
      
      // Check to see if this link density is sufficiently high, compared to the cluster's density.  If so, we have a multi-species contig call.
      if ( link_density_to_cluster < cluster_density[j] * DENSITY_THRESHOLD ) continue;
      
      cout << "Contig #" << i << " is in cluster c" << cID << ", but we think it's in c" << j << " too!\tContig-to-cluster link density: " << link_density_to_cluster << "\tIntra-cluster link density = " << cluster_density[j] << endl;
      
      
      // Determine the reference genomes to which these clusters truly belong.  Then determine whether this contig truly maps to these genomes, which implies
      // whether or not it should be in each of these clusters.  (Note: A single cluster may contain more than one reference genome, e.g., multiple strains of
      // S. cerevisiae.  In this case our choice of one such species is arbitrary, but that's ok because the contig should align to all of the species anyway.)
      int ref0 = _clusters.PluralityRefID( *_truth, cID );
      int ref1 = _clusters.PluralityRefID( *_truth, j );
      if ( ref0 == -1 || ref1 == -1 ) continue; // can't do anything with clusters that have no contigs mapping to any reference
      bool aligns_to_ref0 =_truth->QOnRef( i, ref0 );
      bool aligns_to_ref1 =_truth->QOnRef( i, ref1 );
      
      total_len += _contig_lengths[i];
      
      // The above two Boolean flags imply the "truth" of the multi-species matter.  There are four possibilities:
      // CASE 1: ref0 = true , ref1 = true : True positive!  The contig truly belongs in both species, and thus in both clusters.
      // CASE 2: ref0 = true , ref1 = false: False positive!  The contig only belongs in the species to which it was clustered.
      // CASE 3: ref0 = false, ref1 = true : Clustering error!  The contig shouldn't have been in that other cluster in the first place.  But here we're
      //                                     evaluating the multi-species contig placements, not the clustering, so count this as a true positive.
      // CASE 4: ref0 = false, ref1 = false: No information!  The contig doesn't align anywhere, so who the hell knows what cluster it belongs in.
      if      ( aligns_to_ref1 ) { N_true_positives [cID][j]++; N_true_positives_total++; true_len += _contig_lengths[i]; }
      else if ( aligns_to_ref0 ) { N_false_positives[cID][j]++; N_false_positives_total++; }
      if ( aligns_to_ref1 && !aligns_to_ref0 ) { N_error_fixes_total++; }
      
      
      // Mark Case 1 true positives between the SC and SP clusters.
      if ( aligns_to_ref0 )
	if ( aligns_to_ref1 )
	  if ( ( ref0 == 0 && ref1 == 4 ) || ( ref0 == 4 && ref1 == 0 ) ) {
	    N_refs_04_called++;
	    len_refs_04_called += _contig_lengths[i];
	  }
      
      
    }

    
    
  }
  
  
  
  
  // Find the number and length of multi-species contig calls that SHOULD have been made - that is, all the contigs that might be correctly called as such.
  int N_calls_should = 0;
  int len_calls_should = 0;
  
  
  for ( int i = 0; i < _N_contigs; i++ ) {
    
    // Require contigs to have been clustered.
    if ( _clusters.cluster_ID(i) == -1 ) continue;
    
    // Require contigs to truly belong to multiple species.
    int N_refs_called = 0;
    for ( int j = 0; j < _truth->NRefs(); j++ )
      if ( _truth->QOnRef(i,j) )
	N_refs_called++;
    
    if ( N_refs_called < 2 ) continue;
    
    N_calls_should++;
    len_calls_should += _contig_lengths[i];
  }
  PRINT2( N_calls_should, len_calls_should );
  
  
  
  // Basic report.
  double SENSITIVITY = true_len / double( len_calls_should );
  double SPECIFICITY = N_true_positives_total / double ( N_true_positives_total + N_false_positives_total );
  PRINT7( DENSITY_THRESHOLD, total_len, SENSITIVITY, SPECIFICITY, N_true_positives_total, N_false_positives_total, N_error_fixes_total );
  
  
    
  // Find the number of SC/SP calls that SHOULD have been made.  This is equal to the number of contigs that align to both the SC-FY and SP references, and
  // which were clustered into a cluster whose plurality reference is SC-FY or SP.  These references' IDs are 0 and 4.
  int N_refs_04_should = 0;
  int len_refs_04_should = 0;
  for ( int i = 0; i < _N_contigs; i++ ) {
    if ( !_truth->QOnRef( i, 0 ) ) continue;
    if ( !_truth->QOnRef( i, 4 ) ) continue;
    
    int cID = _clusters.cluster_ID(i);
    if ( cID == -1 ) continue;
    int ref = _clusters.PluralityRefID( *_truth, cID );
    if ( ref != 0 && ref != 4 ) continue;
    
    N_refs_04_should++;
    len_refs_04_should += _contig_lengths[i];
  }
  
  // Find the false negative rate!
  PRINT2(   N_refs_04_called,   N_refs_04_should );
  PRINT2( len_refs_04_called, len_refs_04_should );
}






// Non-truth-based cluster analysis.
// This function uses _link_matrix, _contig_RE_sites, and _clusters.
// Another version of this same calculation is in LinkMatrix::IntraClusterEnrichment(), which inputs the cluster IDs of each contig as a vector.
void
MetaAssembly::FindIntraClusterEnrichment( ostream & out ) const
{
  out << Time() << ": FindIntraClusterEnrichment" << endl;
  
  assert( _link_matrix.has_data() ); // if this fails, you need to load data into the MetaAssembly via LoadFromSAM()
  assert( !_clusters.empty() ); // if this fails, you need to run Cluster() or SetClusters()
  
  int64_t N_intra_cluster_links = 0, N_inter_cluster_links = 0;
  int64_t intra_cluster_len2 = 0, inter_cluster_len2 = 0;
  
  // Loop over all pairs of contigs.  Skip contigs that aren't clustered.
  // This isn't computationally efficient.
  for ( int i = 0; i < _N_contigs; i++ ) {
    int C1 = _clusters.cluster_ID(i);
    if ( C1 == -1 ) continue;
    
    for ( int j = i+1; j < _N_contigs; j++ ) {
      int C2 = _clusters.cluster_ID(j);
      if ( C2 == -1 ) continue;
      
      int N_links = _link_matrix(i,j);
      int len2 = ( _contig_RE_sites[i] - 1 ) * ( _contig_RE_sites[j] - 1 ); // here, remove the +1 that's normally added to contig_RE_sites
      
      ( C1 == C2 ? N_intra_cluster_links : N_inter_cluster_links ) += N_links;
      ( C1 == C2 ?   intra_cluster_len2  :   inter_cluster_len2  ) += len2;
    }
  }
  
  out << Time() << ": Results:" << endl;
  PRINT4( N_intra_cluster_links, N_inter_cluster_links, intra_cluster_len2, inter_cluster_len2 );
  
  double intra_density = double( N_intra_cluster_links ) / double( intra_cluster_len2 );
  double inter_density = double( N_inter_cluster_links ) / double( inter_cluster_len2 );
  
  double enrichment = intra_density / inter_density;
  out << "\t\tIntra-cluster density (links per RE site^2) =\t" << intra_density << endl;
  out << "\t\tInter-cluster density (links per RE site^2) =\t" << inter_density << endl;
  out << "\t\tINTRA-CLUSTER ENRICHMENT =\t" << enrichment << endl;
  
}

 


// Analyze the contigs that aren't clustered.  Gives a more thorough report if truth data is available (_truth != NULL).
void
MetaAssembly::ReportUnclusteredContigs( ostream & out  ) const
{
  out << Time() << ": ReportUnclusteredContigs" << endl;
  
  assert( _link_matrix.has_data() ); // if this fails, you need to load data into the MetaAssembly via LoadFromSAM()
  assert( !_clusters.empty() ); // if this fails, you need to run Cluster() or SetClusters()
  
  // In variable names in this function, "UC" stands for "unclustered contigs", while "CC" stands for "clustered contigs".
  int UC_N = 0;
  vector<int> UC_lengths, UC_RE_sites;
  vector<int> CC_lengths, CC_RE_sites;
  int UC_total_length = 0, UC_total_RE_sites = 0;
  
  // Variables storing reference hits.  Use these only if truth data is available.
  int N_refs = _truth != NULL ? _truth->NRefs() : 0;
  vector<int> UCs_on_ref,    UCs_on_ref_only;
  vector<int> UC_len_on_ref, UC_len_on_ref_only;
  UCs_on_ref        .resize( N_refs, 0 );
  UCs_on_ref_only   .resize( N_refs, 0 );
  UC_len_on_ref     .resize( N_refs, 0 );
  UC_len_on_ref_only.resize( N_refs, 0 );
  int UCs_on_no_ref = 0, UC_len_on_no_ref = 0;
  
  
  // Loop over all contigs and look for UCs (unclustered contigs).
  for ( int i = 0; i < _N_contigs; i++ ) {
    int len = _contig_lengths[i];
    int N_REs = _contig_RE_sites[i] - 1; // here, remove the +1 that's normally added to contig_RE_sites
    
    if ( _clusters.cluster_ID(i) == -1 ) {
      
      // For each UC, keep track of lengths and number of RE sites.
      UC_N++;
      UC_lengths .push_back( len   );
      UC_RE_sites.push_back( N_REs );
      UC_total_length   += len;
      UC_total_RE_sites += N_REs;
      
      // If there's truth data, look at which UCs align (uniquely) to which references.
      if ( _truth != NULL ) {
	for ( int j = 0; j < N_refs; j++ ) {
	  if ( _truth->QOnRef    (i,j) ) { UCs_on_ref     [j]++; UC_len_on_ref     [j] += len; }
	  if ( _truth->QOnRefOnly(i,j) ) { UCs_on_ref_only[j]++; UC_len_on_ref_only[j] += len; }
	}
	
	if ( !_truth->QOnAnyRef(i) ) { UCs_on_no_ref++; UC_len_on_no_ref += len; }
      }
      
      
    }
    
    // Keep track of some statistics for CCs (clustered contigs), for the sake of comparisons to UCs.
    else {
      CC_lengths .push_back( len   );
      CC_RE_sites.push_back( N_REs );
    }
  }
  
  
  // Look for unclustered and unlinked contigs, or "UUCs": contigs that are unclustered AND have no links connecting them to other contigs.
  cout << Time() << ": Looking for contigs that are unclustered and unlinked..." << endl;
  vector<bool> UUC( _N_contigs, true );
  
  for ( int i = 0; i < _N_contigs; i++ ) // mark contigs as clustered...
    if ( _clusters.cluster_ID(i) != -1 )
      UUC[i] = false;
  
  for ( LinkMatrixInt::const_iterator1 it1 = _link_matrix.begin1(); it1 != _link_matrix.end1(); ++it1 ) // ... or linked
    for ( LinkMatrixInt::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2 ) {
      UUC[ it2.index1() ] = false;
      UUC[ it2.index2() ] = false;
    }
  
  // Find contigs that have not been marked as clustered or linked.  Thess are "UUCs".
  int UUC_N = 0, UUC_total_length = 0, UUC_total_RE_sites = 0;
  for ( int i = 0; i < _N_contigs; i++ )
    if ( UUC[i] ) {
      UUC_N++;
      UUC_total_length   += _contig_lengths[i];
      UUC_total_RE_sites += _contig_RE_sites[i] - 1; // here, remove the +1 that's normally added to contig_RE_sites
    }
  
  
  
  // Calculate N50's of all contig length distributions.
  int UC_N50_length    = N50( UC_lengths );
  int UC_N50_RE_sites  = N50( UC_RE_sites );
  int CC_N50_length    = N50( CC_lengths );
  int CC_N50_RE_sites  = N50( CC_RE_sites );
  
  // Calculate percentages.
  double pct_UUC     = 100.0 * UUC_N / UC_N;
  double pct_UUC_len = 100.0 * UUC_total_length   / UC_total_length;
  double pct_UUC_RE  = 100.0 * UUC_total_RE_sites / UC_total_RE_sites;
  
  // Print the report!
  cout << Time() << ": Report!  THERE ARE " << UC_N << " UNCLUSTERED CONTIGS (UC's.)" << endl;
  cout << "UCs' lengths:\tTotal = "  << UC_total_length   << "\tN50 = " << UC_N50_length   << " (vs. " << CC_N50_length   << " for clustered contigs)" << endl;
  cout << "UCs' RE sites:\tTotal = " << UC_total_RE_sites << "\tN50 = " << UC_N50_RE_sites << " (vs. " << CC_N50_RE_sites << " for clustered contigs)" << endl;
  cout << UUC_N << " (" << pct_UUC << "%) of UCs are UUCs, or Unlinked Unclustered Contigs, with no links to other contigs.  Total UUC length = " << UUC_total_length << " (" << pct_UUC_len << "%); total UUC RE sites = " << UUC_total_RE_sites << " (" << pct_UUC_RE << "%)." << endl;
  cout << endl;
  
  if ( _truth != NULL ) {
    cout << "UCs align to the following reference genomes:" << endl;
    for ( int j = 0; j < N_refs; j++ )
      cout << "\tReference #" << j << " = " << _truth->RefName(j) << "\t\tAligning to this: " << UCs_on_ref[j] << " UCs, with length " << UC_len_on_ref[j] << "\t\tAligning ONLY to this: " << UCs_on_ref_only[j] << " UCs, with length " << UC_len_on_ref_only[j] << endl;
    cout << "\tAligning to no reference:\t" << UCs_on_no_ref << " UCs, with length " << UC_len_on_no_ref << endl;
  }
  
}


// DrawFigure2a: Make an igraph output file that depicts this graph and clustering result as an awesome image.
void
MetaAssembly::DrawFigure2a( const string & OUT_PNG ) const
{
  cout << Time() << ": DrawFigure2a" << endl;
  
  assert( _truth != NULL ); // if this fails, you need to call LoadTrueMapping()
  assert( _link_matrix.has_data() ); // if this fails, you need to load data into the MetaAssembly via LoadFromSAM()
  assert( !_clusters.empty() ); // if this fails, you need to run Cluster() or SetClusters()
  
  int N_clusters = _clusters.NClusters();
  
  
  // HEUR: Parameters for graph drawing.
  const string prefix = ""; // "" or "c"
  //const int len_cutoff = 100000;
  const double node_size_scaling_factor = 1e-2; // node size in graph = contig length * node_size_scaling_factor
  const double edge_weight_scaling_factor = 0.5;
  //const int top_N_contigs_per_ref = 100; // a contig will be included in the chart if it's among the top 100 longest contigs that align to its reference...
  //const int top_N_contigs_per_ref_only = 100; // ...or if it's among the top 10 contigs that align only to its reference
  const int top_N_contigs_per_cluster = 10; // TEMP: 200
  const int edge_weight_cutoff = 0; // throw out edges with fewer than this many links
  const double edge_density_cutoff = 0.01; // throw out edges with fewer than this many links per RE_site^2
  const double repulse_exp = 0.2; // 0 = perfect lattice spacing of vertices; 1 = extreme separation of vertices by cluster
  const bool color_by_truth = true; // if 1, only draw contigs that align uniquely to a single species, and color them by true species instead of cluster ID
  
  // Determine which contigs to use as nodes.  Contigs that are too short are not used.
  vector<bool> contigs_to_use( _N_contigs, false );
  
  
  multimap< int,int,greater<int> > contigs_by_len;
  vector< multimap< int,int,greater<int> > > cluster_contigs_by_len( N_clusters );
  
  for ( int i = 0; i < _N_contigs; i++ ) {
    int cID = _clusters.cluster_ID(i);
    contigs_by_len.insert( make_pair( _contig_lengths[i], i ) );
    if ( cID == -1 ) continue; // don't use unclustered contigs
    cluster_contigs_by_len[cID].insert( make_pair( _contig_lengths[i], i ) );
  }
  
  
  
  // This block of code marks the top N contigs aligning to each reference genome for use, regardless of which cluster they're in.
  /*
  vector<int> N_used_on_ref     ( _truth->NRefs(), 0 );
  vector<int> N_used_on_ref_only( _truth->NRefs(), 0 );
  for ( multimap< int,int,greater<int> >::const_iterator it = contigs_by_len.begin(); it != contigs_by_len.end(); ++it ) {
    int contig = it->second;
    
    // If this contig aligns to any reference that has not already seen its allotment of contigs used, then...
    for ( int ref_ID = 0; ref_ID < _truth->NRefs(); ref_ID++ ) {
      if ( !_truth->QOnRef( contig, ref_ID ) ) continue;
      if ( N_used_on_ref     [ref_ID] >= top_N_contigs_per_ref ) continue;
      if ( N_used_on_ref_only[ref_ID] >= top_N_contigs_per_ref_only ) continue;
      
      // ... mark this contig as used, and increment all refs that it aligns to.
      contigs_to_use[contig] = true;
      N_used_on_ref     [ref_ID]++;
      N_used_on_ref_only[ref_ID]++;
    }
    
  }
  */
  
  
  // This block of code marks the top N contigs in each cluster for use.  If color_by_truth = true, it only selects contigs that align uniquely to one ref.
  
  for ( int cID = 0; cID < N_clusters; cID++ ) {
    int N = 0;
    for ( multimap< int,int,greater<int> >::const_iterator it = cluster_contigs_by_len[cID].begin(); it != cluster_contigs_by_len[cID].end(); ++it ) {
      //PRINT4( cID, N, it->first, it->second );
      if ( color_by_truth && _truth->QRefIDOnly( it->second ) < 0 ) continue;
      
      contigs_to_use[ it->second ] = true;
      if ( ++N == top_N_contigs_per_cluster ) break;
    }
  }
  
  
  int N_usable_contigs = count( contigs_to_use.begin(), contigs_to_use.end(), true );
  PRINT( N_usable_contigs );
  
  
  
  
  // Open the igraph_edges file for writing.
  string igraph_edges_file = "Fig2a.edge_chart";
  ofstream igraph_edges_chart( igraph_edges_file.c_str() );
  igraph_edges_chart << "NAME\tNLINKS\tN1\tN2\tCLUSTER\tTRUTH" << endl;
  
  
  
  vector<bool> contigs_used( _N_contigs, false );
  int N_edges = 0;
  
  
  // Loop over all adjacencies in the LinkMatrix.  For each one, write a line to the igraph_edges file describing the edge.
  for ( LinkMatrixInt::const_iterator1 it1 = _link_matrix.begin1(); it1 != _link_matrix.end1(); ++it1 )
    for ( LinkMatrixInt::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2 ) {
      const int c1 = it2.index1(); // contig IDs
      const int c2 = it2.index2();
      
      // TEMP: I'm futzing with these; originally just the first line was decommented
      if ( !contigs_to_use[c1] || !contigs_to_use[c2] ) continue;
      //if ( _contig_lengths[c1] < len_cutoff || _contig_lengths[c2] < len_cutoff ) continue;
      
      // Calculate the edge weight.  If it's not big enough, skip the edge.
      int edge_weight = *it2;
      double edge_density = double(edge_weight) / _contig_RE_sites[c1] / _contig_RE_sites[c2];
      if ( edge_weight  < edge_weight_cutoff )  continue; // if the edge is too small , don't write the edge (cutoff can be total weight, or weight per RE^2)
      if ( edge_density < edge_density_cutoff ) continue;
      
      
      // Write the edge's weight and cluster/truth ID (defined as the same cluster/truth ID as the two nodes, if they're in the same cluster/true ref, and
      // -1 otherwise.)  Note that we're writing the weights in the same order that they appear in the R script.
      int C1 = _clusters.cluster_ID(c1);
      int C2 = _clusters.cluster_ID(c2);
      int cluster = C1 == C2 ? C1 : -1;
      
      int truth1 = _truth->QRefIDOnly(c1);
      int truth2 = _truth->QRefIDOnly(c2);
      int truth = ( truth1 == truth2 && truth1 >= 0 ) ? truth1 : -1;
      
      
      igraph_edges_chart << (N_edges++) << '\t' << edge_weight << "\t" << prefix << c1 << "\t" << prefix << c2 << '\t' << cluster << '\t' << truth << endl;
      
      
      contigs_used[c1] = true;
      contigs_used[c2] = true;
    }
  
  igraph_edges_chart.close();
  
  
  
  // Make a "nodes chart" file containing parameters for each node/contig.
  // The columns are: name; length; cluster ID; truth ID (-1=none, -2=multiple); whether it's used in the plot.
  // These parameters be used by igraph to determine the node/contig's size and color.
  string igraph_nodes_file = "Fig2a.node_chart";
  ofstream igraph_nodes_chart( igraph_nodes_file.c_str(), ios::out );
  igraph_nodes_chart << "NAME\tLENGTH\tCLUSTER\tTRUTH\tUSED" << endl;
  
  for ( int i = 0; i < _N_contigs; i++ )
    igraph_nodes_chart << prefix << i << '\t' << _contig_lengths[i] << '\t' << _clusters.cluster_ID(i) << '\t' << _truth->QRefIDOnly(i) << '\t' << noboolalpha << contigs_used[i] << endl;
  
  igraph_nodes_chart.close();
  
  
  int N_nodes = count( contigs_used.begin(), contigs_used.end(), true );
  PRINT2( N_nodes, N_edges );
  
  // Open the igraph_R file.  This is the R script that we'll be running to create the igraph plot.
  string igraph_R_script = "Fig2a.R";
  cout << Time() << ": Writing R script to file " << igraph_R_script << endl;
  
  ofstream igraph_R( igraph_R_script.c_str(), ios::out );
  
  // Make the beginning of the R script.
  igraph_R <<
    "#!/net/gs/vol3/software/modules-sw/R/2.15.1/Linux/RHEL6/x86_64/bin/Rscript\n"
    "# The above 'shebang' allows this file to be self-executing\n"
    "#\n"
    "# This script was auto-generated by MetaPhase, in the function MetaAssembly::DrawFigure2a().\n"
    "# It uses the data files " << igraph_nodes_file << " and " << igraph_edges_file << ".\n\n\n"
    "library( igraph )\n"
    "library( RColorBrewer )\n\n"
    "print( Sys.time() )\n\n";
  
  
  
  // Write out the part of the R script that reads the igraph_nodes and igraph_edges files (written in this C++ function) and creates tables.
  igraph_R <<
    "node.chart <- read.table( '" << igraph_nodes_file << "', header=TRUE )\n"
    "edge.chart <- read.table( '" << igraph_edges_file << "', header=TRUE )\n\n\n";
  
  igraph_R <<
    "# Make a matrix for the edge list, and use this to load in the graph.\n"
    "edge.list <- cbind( as.character(edge.chart$N1), as.character(edge.chart$N2) )\n"
    "graph <- graph.edgelist( edge.list, directed=FALSE )\n";
  
  
  // Write out the part of the R script that loads in contig lengths and clusters from the igraph nodes chart.
  igraph_R <<
    "# Read in the data about nodes/contigs.  The data has to be remapped by node name.\n"
    "node.names <- V(graph)$name\n"
    "N.nodes <- length( node.names )\n"
    "#node.data <- t( sapply( node.names, function(x) node.chart[node.chart$NAME==x,] ) )\n\n";
  
  // Write out the part of the R script that converts contig clusters into node colors.
  igraph_R <<
    "node.sizes <- sapply( node.names, function(x) node.chart[node.chart$NAME==x,'LENGTH'] )\n"
    "V(graph)$size <- sqrt( node.sizes ) * " << node_size_scaling_factor << "\n\n\n\n";
  
  // Write out the part of the R script that converts contig clusters into node colors.
  string column_name = color_by_truth ? "TRUTH" : "CLUSTER";
  igraph_R <<
    "# Get an array of colors for each node according to the node/contig's cluster.\n"
    "palette <- c( brewer.pal( 3, 'Pastel1' ), brewer.pal( 9, 'Set1' ) )\n" // TODO: de-hardwire number of colors
    "node.colors <- sapply( node.names, function(x) palette[ node.chart[node.chart$NAME==x,'" << column_name << "'] + 1 ] ) # +1 for one-indexing\n"
    "V(graph)$color <- node.colors\n\n\n\n";
  
  // Write out the part of the R script that loads in edge weights and colors from the igraph edges chart.
  // The edge weights and colors will be loaded into the graph object in the plot() function.
  igraph_R <<
    "edge.weights <- sqrt( edge.chart$NLINKS * " << edge_weight_scaling_factor << " )\n\n"
    "# Load edges' cluster IDs and convert them to colors.  This requires a modified palette with an extra color for unclustered edges (-1).\n"
    "palette2 <- c( 'gray20', palette )\n"
    "edge.colors <- palette2[ edge.chart$" << column_name << " + 2 ] # +1 for one-indexing, +1 for the '-1' entries\n"
    "\n\n\n\n";
  
  // Write the rest of the R script that will print out the igraph plot.
  igraph_R <<
    "print( Sys.time() )\n\n"
    "lay <- layout.fruchterman.reingold( graph, niter=500, repulserad=N.nodes*N.nodes*N.nodes^" << repulse_exp << " )\n\n"
    "print( Sys.time() )\n\n"
    "png( file='" << OUT_PNG << "', width=1200, height=1200 )\n\n"
    "plot( graph, layout=lay, vertex.label=NA, edge.width=edge.weights, edge.color=edge.colors )\n"
    "invisible( dev.off() )\n\n"
    "print( Sys.time() )\n\n";
  
  igraph_R.close();
  
  
  // The igraph_R_script file is now a self-executing R script.  Run it!
  cout << Time() << ": Drawing the graph at " << OUT_PNG << endl;
  system( ( "chmod 755 " + igraph_R_script ).c_str() );
  system( igraph_R_script.c_str() );
  
}







// Write this link matrix as an XGMML network file for visualization in Cytoscape.  Nodes are contigs, and edges are Hi-C links between contigs.
// Clustering (if done) and truth information about each contig are included as attributes, and can be used within Cytoscape to change node size, color, etc.
// Contigs are filtered by size (see parameters MIN_N_RES, MIN_LEN, USE_RES); and links between contigs are filtered by frequency (see MIN_N_LINKS).
// The XGMML format is described here: http://wiki.cytoscape.org/XGMML
// Compare to LinkMatrix::WriteXGMML().
void
MetaAssembly::WriteXGMML( const string & XGMML_file, const int MIN_LINK_WEIGHT ) const
{
  cout << Time() << ": WriteXGMML!  Writing to file " << XGMML_file << " for Cytoscape to read..." << endl;
  assert( _link_matrix.has_data() ); // if this fails, you need to load data into the MetaAssembly via LoadFromSAM()
  
  // HEUR: Heuristics for graph appearance.  These are subject to change!
  const int MIN_N_RES = 20;
  const int MIN_LEN = 6000;
  const bool USE_RES = false;
  
  int N_nodes = 0, N_edges = 0;
  
  // Only keep track of fairly large contigs, to keep the graph from getting unwieldy.
  vector<bool> short_contig( _N_contigs, false );
  for ( int i = 0; i < _N_contigs; i++ )
    if ( USE_RES ) short_contig[i] = ( _contig_RE_sites[i] < MIN_N_RES );
    else           short_contig[i] = ( _contig_lengths[i]  < MIN_LEN );
  
  
  // Open the XGMML file and write the XML header.
  ofstream XGMML( XGMML_file.c_str(), ios::out );
  
  XGMML << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<graph label=\"MetaAssembly of awesomeness\" directed=\"0\" xmlns=\"http://www.cs.rpi.edu/XGMML\">" << endl;
  
  
  // For each contig, make a node and write its description.
  for ( int i = 0; i < _N_contigs; i++ ) {
    if ( short_contig[i] ) continue;
    
    N_nodes++;
    
    XGMML << "  <node id=\"" << i << "\" label=\"" << _contig_names[i] << "\">" << endl;
    XGMML << "    <att name=\"size\" type=\"integer\" value=\"" << _contig_lengths[i] << "\"/>" << endl;
    
    // Describe the contig's TrueMapping, if known.
    if ( _truth != NULL ) {
      XGMML << "    <att name=\"refs_bitcode\" type=\"integer\" value=\"" << _truth->RefsBitcode(i) << "\"/>" << endl;
      // The TrueMapping includes two flags for each reference: one if the contig aligns to that ref, and one if the contig aligns ONLY to that ref.
      for ( int j = 0; j < _truth->NRefs(); j++ ) {
	XGMML << "    <att name=\"on_ref_"      << _truth->RefName(j) << "\" type=\"integer\" value=\"" << _truth->QOnRef    (i,j) << "\"/>" << endl;
	XGMML << "    <att name=\"only_on_ref_" << _truth->RefName(j) << "\" type=\"integer\" value=\"" << _truth->QOnRefOnly(i,j) << "\"/>" << endl;
      }
    }
    
    // Describe the contig's clustering result, if known.
    if ( !_clusters.empty() )
      XGMML << "    <att name=\"cluster\" type=\"integer\" value=\"" << _clusters.cluster_ID(i) << "\"/>" << endl;
    
    XGMML << "  </node>\n";
  }
  
  // For each adjacency, write an edge.
  
  for ( int i = 0; i < _N_contigs; i++ ) {
    if ( short_contig[i] ) continue;
    
    for ( int j = i+1; j < _N_contigs; j++ ) { // note that j > i
      if ( short_contig[j] ) continue;
      
      double link_weight = _link_matrix(i,j);
      //link_weight /= ( _contig_RE_sites[i] * _contig_RE_sites[j] );
      if ( link_weight < MIN_LINK_WEIGHT ) continue;
      
      N_edges++;
      
      XGMML << "  <edge source=\"" << i << "\" target=\"" << j << "\">" << endl;
      XGMML << "    <att name=\"link_weight\" type=\"real\" value=\"" << link_weight << "\"/>" << endl;
      XGMML << "  </edge>" << endl;
    }
  }
  
  
  XGMML << "</graph>" << endl;
  XGMML.close();
  
  
  // Write a summary of what just happened.
  cout << Time() << ": WriteXGMML: Done!  Wrote " << N_nodes << " nodes (representing contigs " << ( USE_RES ? "with N RE sites >= " : "of length >= " )
       << ( USE_RES ? MIN_N_RES : MIN_LEN ) << ") and " << N_edges << " edges (each representing >=" << MIN_LINK_WEIGHT << " Hi-C links between contigs)."
       << endl;
}



// Write a fasta file representing all the contigs in one cluster.  The cluster chosen is cluster #ID - or, if is_ref_ID = true, the single largest cluster
// corresponding to ref ID #ID.  If ID == -1 and !is_ref_ID, writes all of the contigs not placed in any cluster.
void
MetaAssembly::WriteClusterFasta( const string & output_fasta, const int ID, const bool is_ref_ID ) const
{
  cout << Time() << ": WriteClusterFasta: " << ( is_ref_ID ? "cluster for ref #" : "cluster #" ) << ID << "\t->\t" << output_fasta << endl;
  
  assert( !_clusters.empty() );
  assert( boost::filesystem::is_regular_file( _assembly_fasta ) ); // technically this isn't quite required by the constructor, so we're double-checking
  
  
  
  // Find the cluster ID that we will be using.  If is_ref_ID, we need to find the cluster ID that corresponds to this reference ID.
  int cID;
  
  if ( !is_ref_ID ) {
    assert( ID >= -1 ); // if !is_ref_ID, ID == -1 is allowable, and indicates unclustered contigs
    assert( ID < _clusters.NClusters() );
    
    cID = ID;
  }
  
  else {
    assert( ID >= 0 );
    assert( ID < _truth->NRefs() );
    
    // Find the largest cluster whose plurality reference is this reference ID.  This is the cluster that we will output.
    cID = _clusters.LargestClusterWithPluralityRefID( *_truth, ID );
    
    // If no such cluster could be found, there's nothing to do.
    if ( cID == -1 ) {
      cout << "Oops, we couldn't find a single cluster whose plurality reference was " << ID << " = " << _truth->RefName(ID) << ".  So we can't write any fasta files. Sorry!" << endl;
      return;
    }
  }
  
  
  // Create the directory to contain this output file, if necessary.
  int pos = output_fasta.rfind("/");
  string output_dir = output_fasta.substr( 0, pos );
  system( ( "mkdir -p " + output_dir ).c_str() );
  
  
  // Method: Read the assembly fasta file, and write to the the cluster fasta file.
  static const unsigned LINE_LEN = 1000000;
  char LINE[LINE_LEN];
  
  ifstream in( _assembly_fasta.c_str(), ios::in );
  ofstream out(   output_fasta.c_str(), ios::out);
  
  int contig_ID = -1; // contig's order in assembly fasta; also same as the contig ID in the ClusteringResult; starts at -1 because it will be incremented to 0
  
 
  // Read the assembly fasta file line-by-line.
  while ( 1 ) {
    in.getline( LINE, LINE_LEN );
    assert( strlen(LINE) + 1 < LINE_LEN );
    if ( in.fail() ) break;
    
    // Check to see if this line starts a new contig.  In this way, keep track of the contig ID corresponding to every line in the fasta file.
    if ( LINE[0] == '>' ) contig_ID++;
    
    // If this contig belongs in this cluster (including the possibility of -1 indicating no clusters), write it to the output file.
    if ( _clusters.cluster_ID( contig_ID ) == cID ) out << LINE << endl;
  }
  
  assert( contig_ID + 1 == _N_contigs ); // if this fails, the assembly fasta might have been parsed wrong or might be the wrong fasta
  
  
  in.close();
  out.close();
}




// Write a set of files describing one cluster, to be input to Lachesis.  The cluster chosen is the single largest cluster corresponding to ref #ref_ID.
// The files include: all.GLM only (that should be all that's needed; in the future maybe pare down the SAM file?)
void
MetaAssembly::WriteLachesisFiles( const string & Lachesis_dir, const int ref_ID ) const
{
  cout << Time() << ": WriteLachesisFiles ->\t" << Lachesis_dir << endl;
  
  assert( !_clusters.empty() );
  assert( ref_ID >= 0 );
  assert( ref_ID < _truth->NRefs() );
  
  // Find the largest cluster whose plurality reference is this reference ID.  This is the cluster that we will output.
  int cID = _clusters.LargestClusterWithPluralityRefID( *_truth, ref_ID );
  
  // If no such cluster could be found, there's nothing to do.
  if ( cID == -1 ) {
    cout << "Oops, we couldn't find a single cluster whose plurality reference was " << ref_ID << " = " << _truth->RefName(ref_ID) << ".  So we can't write any Lachesis files. Sorry!" << endl;
    return;
  }
  
  
  // Set up the output directory.
  system( ( "mkdir -p " + Lachesis_dir ).c_str() );
  
  string GLM_file = Lachesis_dir + "/all.GLM";
  ofstream out( GLM_file.c_str(), ios::out );
  
  // Write the all.GLM file, which describes the non-normalized matrix of Hi-C links.
  // This file starts with a header, whose format is prescribed by Lachesis.  See the functions ReadFile() and WriteFile() in Lachesis's GenomeLinkMatrix.cc.
  out << "# GenomeLinkMatrix file - produced by MetaPhase's MetaAssembly::WriteLachesisFiles() - see GenomeLinkMatrix.h for documentation of this object type" << endl;
  out << "# Species = " << _truth->RefName(ref_ID) << endl;
  out << "# N_bins = " << _N_contigs << endl;
  out << "# bin_size = 0" << endl;
  out << "# RE_sites_file = " << _RE_sites_file << endl;
  out << "# SAM files used in generating this dataset:";
  for ( size_t i = 0; i < _SAM_files.size(); i++ )
    out << " " << _SAM_files[i];
  out << endl;
  
  
  int N_links_written = 0;
  
  // Print the non-normalized link counts to file.  However, ONLY print links between contigs that have been clustered into our cluster of interest.
  // The idea is that alhough Lachesis will see all of the contigs in the entire metagenome draft assembly, it will only see the linkages in this MetaPhase
  // cluster, and thus will only create chromosome clusters and orderings for the contigs that MetaPhase has determined to be in this particular species.
  for ( LinkMatrixInt::const_iterator1 it1 = _link_matrix.begin1(); it1 != _link_matrix.end1(); ++it1 )
    for ( LinkMatrixInt::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2 ) {
      if ( _clusters.cluster_ID( it2.index1() ) != cID ) continue;
      if ( _clusters.cluster_ID( it2.index2() ) != cID ) continue;
      
      // Note that MetaPhase stores its LinkMatrixInt as an upper triangular matrix, but Lachesis expects a symmetric matrix.
      out << it2.index1() << '\t' << it2.index2() << '\t' << *it2 << endl;
      out << it2.index2() << '\t' << it2.index1() << '\t' << *it2 << endl;
      N_links_written++;
    }
  
  out.close();
  
  cout << Time() << ": Wrote the " << N_links_written << " links internal to cluster #" << cID << "." << endl;
  
}




