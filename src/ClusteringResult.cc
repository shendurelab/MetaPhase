// For documentation, see ClusteringResult.h
#include "ClusteringResult.h"
#include "TrueMapping.h"
#include "HierarchicalClustering.h" // AgglomerativeHierarchicalClusteringSimple
#include "LinkMatrix.h" // LinkMatrixDouble
#include "TextFileParsers.h" // TokenizeFile

#include <assert.h>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <iostream>
#include <numeric> // accumulate
#include <algorithm> // count, max_element, distance

// Boost includes
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/lexical_cast.hpp>

// Modules in ~/include (must add -L~/include and -lJ<module> to link)
#include "TimeMem.h"





// Null constructor.
ClusteringResult::ClusteringResult()
{
  _N_contigs = 0;
  _N_clusters = 0;
}



// Construct a ClusteringResult with a set of contigs but no cluster assignments.
// Allocate space for clusters if desired, though.
ClusteringResult::ClusteringResult( const vector<int> & contig_lengths, const int N_clusters )
{
  _N_contigs = contig_lengths.size();
  _contig_lengths = contig_lengths;
  _contig_clusters = vector<int>( _N_contigs, -1 );
  
  _N_clusters = N_clusters;
  _cluster_lengths = vector<int>( N_clusters, 0 );
}



// Construct a ClusteringResult object from the immediate result of clustering.
ClusteringResult::ClusteringResult( const vector<int> & contig_to_clusterID, const vector<int> & contig_lengths )
{
  assert( _contig_clusters.size() == _contig_lengths.size() );
  _contig_clusters = contig_to_clusterID;
  _contig_lengths = contig_lengths;
  
  _N_contigs = _contig_clusters.size();
  assert( _N_contigs > 0 );
  
  
  // Find the number of clusters from the maximum value of _contig_clusters.  Then find cluster lengths in bp.
  // NOTE: this assumes the clusters aren't skip-numbered!
  _N_clusters = *( max_element( _contig_clusters.begin(), _contig_clusters.end() ) ) + 1;
  FindClusterLengths();
}







// ReadFile, WriteFile: Read/write a clustering result from file.
// This alows you to evade the runtime cost of running the MetaAssembly::Cluster() function, once you've already run it with the same dataset and parameters.
void
ClusteringResult::ReadFile( const string & filename )
{
  // Load the entire file into tokens.
  vector< vector<string> > tokens;
  TokenizeFile( filename, tokens, true );
  
  // Determine the number of contigs.
  _N_contigs = tokens.size();
  assert( _N_contigs > 0 );
  _contig_lengths .resize( _N_contigs, -1 );
  _contig_clusters.resize( _N_contigs, -1 );
  
  // Each line in the input file should have two tab-delimited tokens: the contig length and the contig's cluster ID.
  for ( int i = 0; i < _N_contigs; i++ ) {
    assert( tokens[i].size() == 2 );
    _contig_lengths [i] = boost::lexical_cast<int>( tokens[i][0] );
    _contig_clusters[i] = boost::lexical_cast<int>( tokens[i][1] );
  }
  
  // Find the number of clusters from the maximum value of _contig_clusters.  Then find cluster lengths in bp.
  _N_clusters = *( max_element( _contig_clusters.begin(), _contig_clusters.end() ) ) + 1;
  FindClusterLengths();
}



void
ClusteringResult::WriteFile( const string & filename ) const
{
  assert( !_contig_clusters.empty() ); // if this fails, you haven't yet run Cluster()
  assert( _N_contigs == (int) _contig_clusters.size() );
  
  ofstream out( filename.c_str() );
  for ( int i = 0; i < _N_contigs; i++ )
    out << _contig_lengths[i] << '\t' << _contig_clusters[i] << endl;
  out.close();
}






// Add a contig to a cluster.  Remove it from whatever cluster it was in before.  This function can be used to create a new cluster at #_N_clusters.
void
ClusteringResult::SetCluster( const int contig_ID, const int cluster_ID )
{
  assert( contig_ID >= 0 && contig_ID < _N_contigs );
  assert( cluster_ID >= 0 && cluster_ID <= _N_clusters );
  
  // Create a new cluster, if necessary.
  if ( cluster_ID == _N_clusters ) {
    _N_clusters++;
    _cluster_lengths.push_back(0);
  }
  
  // If this contig was already assigned to a cluster, remove it from that cluster.
  if ( _contig_clusters[contig_ID] != -1 ) _cluster_lengths[ _contig_clusters[contig_ID] ] -= _contig_lengths[contig_ID];
  
  // Assign the contig to the new cluster.
  _contig_clusters[contig_ID] = cluster_ID;
  _cluster_lengths[cluster_ID] += _contig_lengths[contig_ID];
}
  





// ReorderClustersByRefs: Renumber the cluster IDs (implied by _contig_clusters and _cluster_lengths) in increasing order of PluralityRefID().
// This uses truth information, but the cluster order doesn't really matter and is just for convenience/visual appeal, so this isn't really cheatin
void
ClusteringResult::ReorderClustersByRefs( const TrueMapping & truth )
{
  cout << Time() << ": ClusteringResult::ReorderClustersByRefs" << endl;
  if ( _N_contigs != (int) truth.NQuerySeqs() ) PRINT2( _N_contigs, truth.NQuerySeqs() );
  assert( _N_contigs == (int) truth.NQuerySeqs() );
  
  
  
  // Find the "plurality reference" for each cluster - that is, the reference genome ID to which a plurality of that cluster's sequence aligns.
  // Put these numbers into a map to sort them by plurality ref ID.
  multimap<int,int> plurality_ref_to_cID;
  for ( int i = 0; i < _N_clusters; i++ )
    plurality_ref_to_cID.insert( make_pair( PluralityRefID( truth,i ), i ) );
  
  // Use the map to determine the new ordering of cluster IDs.
  vector<int> old_to_new_cID( _N_clusters, -1 );
  int i = 0;
  for ( multimap<int,int>::const_iterator it = plurality_ref_to_cID.begin(); it != plurality_ref_to_cID.end(); ++it )
    //PRINT2( it->second, it->first );
    old_to_new_cID[ it->second ] = i++;
  
  
  
  // Apply the new ordering to _contig_clusters and _cluster_lengths.
  assert( (int) old_to_new_cID.size() == _N_clusters );
  for ( int i = 0; i < _N_contigs; i++ )
    if ( _contig_clusters[i] != -1 )
      _contig_clusters[i] = old_to_new_cID[ _contig_clusters[i] ];
  
  vector<int> new_cluster_lengths( _N_clusters );
  for ( int i = 0; i < _N_clusters; i++ )
    new_cluster_lengths[ old_to_new_cID[i] ] = _cluster_lengths[i];
  _cluster_lengths = new_cluster_lengths;
}



// ClusterSize: return the number of contigs in this cluster (currently computationally inefficient for multiple calls)
int
ClusteringResult::ClusterSize( const int cluster_ID ) const
{
  assert( cluster_ID >= 0 && cluster_ID < _N_clusters );
  return count( _contig_clusters.begin(), _contig_clusters.end(), cluster_ID );
}




// PluralityRefID: Determine the reference ID to which a plurality of the contigs in this cluster align.  This uses truth info so it's a kind of cheating.
// Returns -1 if no contig in this cluster aligns to any reference.
int
ClusteringResult::PluralityRefID( const TrueMapping & truth, const int cluster_ID ) const
{
  assert( cluster_ID >= 0 );
  assert( cluster_ID < _N_clusters );
  assert( _N_contigs == (int) truth.NQuerySeqs() );
  
  
  // Determine the amount of sequence in this cluster that aligns to each reference genome.
  int N_refs = truth.NRefs();
  vector<int> len_on_ref( N_refs, 0 );
  
  for ( int i = 0; i < _N_contigs; i++ ) {
    if ( _contig_clusters[i] != cluster_ID ) continue;
    
    for ( int j = 0; j < N_refs; j++ )
      if ( truth.QOnRef(i,j) )
	len_on_ref[j] += _contig_lengths[i];
  }
  
  // Whichever reference has the most sequence from this cluster is the plurality reference.
  int plurality_ref = max_element( len_on_ref.begin(), len_on_ref.end() ) - len_on_ref.begin();
  if ( plurality_ref == 0 && len_on_ref[0] == 0 ) return -1;
  return plurality_ref;
}
  

// PluralityRefName: Wrapper for PluralityRefID.  Returns "" if no contig in this cluster aligns to any reference.
string
ClusteringResult::PluralityRefName( const TrueMapping & truth, const int cluster_ID ) const
{
  int plurality_ref_ID = PluralityRefID( truth, cluster_ID );
  if ( plurality_ref_ID == -1 ) return "";
  return truth.RefName( plurality_ref_ID );
}



// Find the cluster ID of the largest cluster whose plurality reference is this reference ID.  "Largest cluster" is measured by contig length.
// Return -1 if no cluster has this plurality ref ID.
int
ClusteringResult::LargestClusterWithPluralityRefID( const TrueMapping & truth, const int ref_ID ) const
{
  assert( ref_ID >= 0 && ref_ID < truth.NRefs() );
  
  int cID = -1;
  int cID_length = 0;
  for ( int i = 0; i < _N_clusters; i++ ) {
    if ( PluralityRefID( truth, i ) != ref_ID ) continue; // TODO: this could be more computationally efficient
    
    if ( _cluster_lengths[i] > cID_length ) {
      cID_length = _cluster_lengths[i];
      cID = i;
    }
  }
  
  return cID; // might still be -1
}









// Print: Output summary statistics.
void
ClusteringResult::Print( const TrueMapping & truth, ostream & out ) const
{
  out << Time() << ": ClusteringResult::Print!" << endl;
  
  PRINT( _N_clusters );
  
  // Variables that will hold statistics.
  vector<int> N_contigs_per_cluster( _N_clusters, 0 );
  int N_clustered_contigs = 0;
  int len_clustered_contigs = 0;
  int N_unclustered_contigs = 0;
  int len_unclustered_contigs = 0;
  int N_misclustered_contigs = 0;
  int len_misclustered_contigs = 0;
  
  
  // Loop over all contigs and collect statistics.
  for ( int i = 0; i < _N_contigs; i++ ) {
    int cluster = _contig_clusters[i];
    
    if ( cluster == -1 ) { // unclustered contigs
      N_unclustered_contigs++;
      len_unclustered_contigs += _contig_lengths[i];
      continue;
    }
    
    // Clustered contigs.
    N_clustered_contigs++;
    len_clustered_contigs += _contig_lengths[i];
    N_contigs_per_cluster[cluster]++;
    int cluster_ref_ID = PluralityRefID( truth, cluster );
    if ( cluster_ref_ID == -1 ) continue;
    
    // The criterion for a mis-clustered contig is: (1) it maps to at least one reference, and (2) it has been clustered, but (3) it doesn't map to the
    // plurality reference in its cluster.
    if ( truth.RefsBitcode(i) != 0 && !truth.QOnRef( i, cluster_ref_ID ) ) {
      //PRINT4( i, cluster, cluster_ref_ID, truth.RefsBitcode(i) );
      N_misclustered_contigs++;
      len_misclustered_contigs += _contig_lengths[i];
    }
  }
  
  assert( _N_contigs == N_clustered_contigs + N_unclustered_contigs );
  
  // Calculate percentages.
  double pct_clustered        = 100.0 * double(     N_clustered_contigs) / _N_contigs;
  double pct_len_clustered    = 100.0 * double(   len_clustered_contigs) / ( len_clustered_contigs + len_unclustered_contigs );
  double pct_misclustered     = 100.0 * double(  N_misclustered_contigs) /   N_clustered_contigs;
  double pct_len_misclustered = 100.0 * double(len_misclustered_contigs) / len_clustered_contigs;
  
  
  
  
  // Report!
  for ( int i = 0; i < _N_clusters; i++ ) {
    cout << "Cluster #" << i << " has " << N_contigs_per_cluster[i] << "\tcontigs, with total length " << _cluster_lengths[i] << ".\tIts plurality reference is #" << PluralityRefID(truth,i) << " (" << PluralityRefName(truth,i) << ")." << endl;
  }
  
  
  cout << "RESULT: " << N_clustered_contigs   << " contigs (" << pct_clustered << "%) are assigned to clusters; total length = " << len_clustered_contigs << " (" << pct_len_clustered << "%).\t" << N_misclustered_contigs << " (" << pct_misclustered << "%) contigs are misclustered; total length = " << len_misclustered_contigs << " (" << pct_len_misclustered << "%)." << endl;
  cout << "(Misclustered contigs are defined as contigs that are clustered, and map to at least one reference, but not to their cluster's plurality reference.)" << endl;
  
  cout << "There are " << N_unclustered_contigs << " contigs not assigned to any cluster.  Total length = " << len_unclustered_contigs << endl;
} 






// Truth-based evaluation of the clusters in this ClusteringResult.  Produce the truth heatmap and chart files.  If plot = true, plot the heatmap image too.
// If require_unique_alignment = true, only count contigs aligning uniquely to a reference (single-species contigs)
// Only draw clusters with length (in bp) at least MIN_CLUSTER_LEN (default 0).
void
ClusteringResult::DrawChart( const TrueMapping & truth, const bool require_unique_alignment, const bool plot, const bool MY, const int MIN_CLUSTER_LEN ) const
{
  cout << Time() << ": ClusteringReport::DrawChart" << endl;
  
  assert( _N_contigs == truth.NQuerySeqs() );
  
  
  // Options
  bool MERGE_SC = MY; // skip SC-non-FY, condense the SC strains
  bool HIDE_DEAD = MY;
  if ( MERGE_SC || HIDE_DEAD ) assert( truth.NRefs() == 16 && truth.RefName(1) == "SC-CEN" ); // make sure this is the MetaYeast scenario!
  
  
  // For each cluster, and each reference genome, count the number of contigs and the amount of sequence in that cluster that derives from that genome (only).
  int N_refs = truth.NRefs();
  const vector<int> init( N_refs, 0 );
  vector< vector<int> > cluster_N_to_ref       ( _N_clusters, init );
  vector< vector<int> > cluster_len_to_ref     ( _N_clusters, init );
  vector< vector<int> > cluster_N_to_ref_only  ( _N_clusters, init );
  vector< vector<int> > cluster_len_to_ref_only( _N_clusters, init );
  
  // Loop over all contigs that are on sufficiently large clusters.
  for ( int i = 0; i < _N_contigs; i++ ) {
    int cluster = _contig_clusters[i];
    if ( cluster == -1 ) continue;
    
    int len = _contig_lengths[i];
    
    // Loop over all references.
    for ( int j = 0; j < N_refs; j++ ) {
      if ( truth.QOnRef(i,j) ) {
	cluster_N_to_ref  [cluster][j]++;
	cluster_len_to_ref[cluster][j] += len;
      }
      if ( MERGE_SC ) { // Consider contigs that align only to the S. cerevisiae strains as "unique alignments"
	bool only = false;
	if ( j >= 4 && truth.QOnRefOnly(i,j) ) only = true;
	if ( j <  4 && truth.RefsBitcode(i) <= 15 ) only = true;
	if ( only ) {
	  cluster_N_to_ref_only  [cluster][j]++;
	  cluster_len_to_ref_only[cluster][j] += len;
	}
      }
      else
      if ( truth.QOnRefOnly(i,j) ) {
	cluster_N_to_ref_only  [cluster][j]++;
	cluster_len_to_ref_only[cluster][j] += len;
      }
    }
  }
  
  
  
  // Write the results to a heatmap file (out), and also to a text chart file (out_chart).
  string suffix = "ClusteringResult."
    + ( require_unique_alignment ? string("u.") : string("nu.") )
    + boost::lexical_cast<string>( NClusters() ) + ".txt";
  ofstream out_chart  ( ( "chart."   + suffix ).c_str(), ios::out ); 
  ofstream out_heatmap( ( "heatmap." + suffix ).c_str(), ios::out );
  
  
  // Header line for chart file.
  out_chart << "\t\t      LENGTH";
  for ( int j = 0; j < N_refs; j++ ) {
    if ( MERGE_SC )  if ( j == 1 || j == 2 || j == 3 ) continue;
    if ( HIDE_DEAD ) if ( j == 15 ) continue;
    out_chart << '\t' << truth.RefName(j);
  }
  out_chart << endl;
  
  
  int y = 0; // y-value of points in heatmap
    
  for ( int i = 0; i < _N_clusters; i++ ) {
    if ( _cluster_lengths[i] < MIN_CLUSTER_LEN ) continue; // save space by not printing lines for too-small clusters
    y++;
    
    char x[12];
    sprintf( x, "%12d", _cluster_lengths[i] );
    out_chart << "CLUSTER #" << i << "\t" << x;
    
    for ( int j = 0; j < N_refs; j++ ) {
      if ( MERGE_SC )  if ( j == 1 || j == 2 || j == 3 ) continue;
      if ( HIDE_DEAD ) if ( j == 15 ) continue;
      int len = require_unique_alignment ? cluster_len_to_ref_only[i][j] : cluster_len_to_ref[i][j];
      double pct_to_ref = 100.0 * double(len) / _cluster_lengths[i];
      out_heatmap << ( MERGE_SC && j == 0 ? "S.cerevisiae" : truth.RefFullName(j) ) << '\t' << y << '\t' << pct_to_ref << endl;
      sprintf( x, "%.2f", pct_to_ref );
      out_chart << '\t' << x << "%";
    }
    out_chart << endl;
  }
  
  
  out_chart  .close();
  out_heatmap.close();
  
  
  if ( plot ) system ( ( "QuickHeatmap heatmap." + suffix ).c_str() );
}






 

// DrawFigure2bc: Call DrawChart to create the chart and heatmap files, then call a script in ../figs to generate Figures 2b and 2c of the paper.
void
ClusteringResult::DrawFigure2bc( const TrueMapping & truth, const bool MY ) const
{
  // First call DrawChart to create the files heatmap.ClusteringResult.A.B.txt, where A = "u", "nu" (unique, nonunique) and B = NClusters().
  DrawChart( truth, false, false, MY, 500000 );
  DrawChart( truth, true,  false, MY, 500000 );
  
  // Now call the script that makes the heatmaps.
  for ( int i = 0; i < 2; i++ ) {
    string version = i ? "nu" : "u";
    string fignum = i ? "2c" : "2b";
    string file = "heatmap.ClusteringResult." + version + "." + boost::lexical_cast<string>( NClusters() ) + ".txt";
    system ( ( "../figs/2/MakeClusteringResultHeatmap.R " + file ).c_str() );
    system ( ( "cp ~/public_html/MakeClusteringResultHeatmap.jpg ~/public_html/Fig" + fignum + ".jpg" ).c_str() ); // put the heatmaps in place
  }
  
} 




// FindClusterLengths: Input _contig_lengths and _contig_clusters; fill _cluster_lengths.  This is purely a bookkeeping operation.
void
ClusteringResult::FindClusterLengths()
{
  _cluster_lengths.resize( _N_clusters );
  
  // Find the sum length of all contigs in each cluster.
  for ( int i = 0; i < _N_contigs; i++ ) {
    int cluster = _contig_clusters[i];
    if ( cluster == -1 ) continue;
    //cluster_contigs[cluster].insert( i ); 
    _cluster_lengths[cluster] += _contig_lengths[i];
  }
}




/* MergeClusteringResultsX: Combine two sets of ClusteringResults into a single ClusteringResult by matching up the clusters and finding intelligent ways to
 * split/merge/hem/haw little bits of cluster.
 * These functions are designed as an intelligent way to merge ClusteringResults produced by Hi-C datasets with differing restriction sites.
 *
 * Method (MergeClusteringResults1): straightforward "intersection" method.
 * 1. Find "matches" between clusters.  That is, for each cluster in each ClusteringResult, find which cluster(s) in the other ClusteringResult seem to
 *    contain the same points.  These matches form a bipartite graph.
 * 2. Create a Boost graph object for this bipartite graph, and find connected components.
 * 3. For each pair of clusters that match uniquely to each other (form a component), prune their clusters in the result.
 */
ClusteringResult
MergeClusteringResults1( const ClusteringResult & cr1, const ClusteringResult & cr2 )
{
  cout << Time() << ": MergeClusteringResults1!" << endl;
  
  static const int MIN_CLUSTER_LEN = 100000;
  
  // Make sure these ClusteringResults are on the same dataset!
  assert( cr1._N_contigs      == cr2._N_contigs );
  assert( cr1._contig_lengths == cr2._contig_lengths );
  
  const int N_contigs = cr1._N_contigs;
  const int N1 = cr1._N_clusters, N2 = cr2._N_clusters;
  
  // Find the sizes of each cluster in bp.  We will ignore all clusters of size < MIN_CLUSTER_LEN.
  vector<int> cr1_sizes, cr2_sizes;
  for ( int c1 = 0; c1 < N1; c1++ )
    cr1_sizes.push_back( cr1.ClusterLength(c1) );
  for ( int c2 = 0; c2 < N2; c2++ )
    cr2_sizes.push_back( cr2.ClusterLength(c2) );
  
  
  
  // 1. Find "matches" between clusters.  That is, for each cluster in each ClusteringResult, find which cluster(s) in the other ClusteringResult seem to
  //    contain the same points.  These matches form a bipartite graph.
  
  // First, determine the amount of sequence length that is in cluster #c1 in cr1 and in cluster #c2 in cr2, for all values of i,j.
  vector< vector<int> > seq_in_c1c2( N1, vector<int>( N2, 0 ) );
  vector< vector<int> > seq_in_c2c1( N2, vector<int>( N1, 0 ) );
  for ( int i = 0; i < N_contigs; i++ ) {
    int c1 = cr1._contig_clusters[i];
    int c2 = cr2._contig_clusters[i];
    if ( c1 == -1 || c2 == -1 ) continue; // ignore contigs that aren't placed in both clusters
    //if ( cr1_sizes[c1] < MIN_CLUSTER_LEN || cr2_sizes[c2] < MIN_CLUSTER_LEN ) continue;
    
    seq_in_c1c2[c1][c2] += cr1._contig_lengths[i];
    seq_in_c2c1[c2][c1] += cr1._contig_lengths[i];
  }
  
  // Now, for each cluster in both ClusteringResults, find the corresponding cluster in the other ClusteringResult with the plurality of its sequence.
  vector<int> match_c1_to_c2( N1 );
  vector<int> match_c2_to_c1( N2 );
  vector< set<int> > match_c1_from_c2( N2 );
  vector< set<int> > match_c2_from_c1( N1 );
  
  for ( int c1 = 0; c1 < N1; c1++ ) {
    int c2 = distance( seq_in_c1c2[c1].begin(), ( max_element( seq_in_c1c2[c1].begin(), seq_in_c1c2[c1].end() ) ) );
    match_c1_to_c2[c1] = c2;
    match_c2_from_c1[c2].insert( c1 );
    cout << "cr1:" << c1 << "\t->\tcr2:" << c2 << endl;
  }
  for ( int c2 = 0; c2 < N2; c2++ ) {
    int c1 = distance( seq_in_c2c1[c2].begin(), ( max_element( seq_in_c2c1[c2].begin(), seq_in_c2c1[c2].end() ) ) );
    match_c2_to_c1[c2] = c1;
    match_c1_from_c2[c1].insert( c2 );
    cout << "cr1:" << c1 << "\t<-\tcr2:" << c2 << endl;
  }
  
  
  // Alternative method.  Makes more clusters, but potentially better ones.
  map< pair<int,int>, int > c1c2_to_new_cID;
  int new_cID = 0;
  for ( int c1 = 0; c1 < N1; c1++ )
    for ( int c2 = 0; c2 < N2; c2++ )
      if ( seq_in_c1c2[c1][c2] >= MIN_CLUSTER_LEN ) {
	PRINT3( c1, c2, seq_in_c1c2[c1][c2] );
	c1c2_to_new_cID.insert( make_pair( make_pair(c1,c2), new_cID++ ) );
      }
  
  ClusteringResult merged_clusters( cr1._contig_lengths, new_cID );
  
  for ( int i = 0; i < N_contigs; i++ ) {
    int c1 = cr1._contig_clusters[i];
    int c2 = cr2._contig_clusters[i];
    if ( c1 == -1 || c2 == -1 ) continue; // ignore contigs that aren't placed in both clusters
    
    if ( seq_in_c1c2[c1][c2] < MIN_CLUSTER_LEN ) continue; // ignore cluster pairings that don't have a lot of sequence in them
    int new_cID = c1c2_to_new_cID[ make_pair(c1,c2) ];
    //PRINT4( i, c1, c2, new_cID );
    merged_clusters.SetCluster( i, new_cID );
  }
  
  return merged_clusters;
  
  
  
  // 2. Create a Boost graph object for this bipartite graph, and find connected components.
  // Note that the graph is undirected even though the links can be directed.
  typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS, boost::no_property > Graph;
  Graph G( N1 + N2 );
  for ( int i = 0; i < N1; i++ )
    add_edge( i, N1 + match_c1_to_c2[i], G );
  for ( int i = 0; i < N2; i++ )
    add_edge( N1 + i, match_c2_to_c1[i], G );
  
  vector<int> CCs( N1 + N2 );
  int N_components = connected_components( G, &CCs[0] );
  N_components++;
  
  for ( int i = 0; i < N1 + N2; i++ )
    PRINT2( i, CCs[i] );
  
  
  
  
  // 3. For each pair of clusters that match uniquely to each other (form a component), prune their clusters in the result.
  
  for ( int c1 = 0; c1 < N1; c1++ ) {
    int c2 = match_c1_to_c2[c1];
    if ( match_c1_from_c2[c1].size() != 1 || match_c2_from_c1[c2].size() != 1 ) continue;
    
    PRINT2( c1, c2 );
    
    for ( int i = 0; i < N_contigs; i++ )
      if ( merged_clusters._contig_clusters[i] != c1 && cr2._contig_clusters[i] == c2 ) {
	cout << "REMOVING: " << i << endl;
	merged_clusters._contig_clusters[i] = -1;
	merged_clusters._cluster_lengths[c1] -= merged_clusters._contig_lengths[i];
      }
  }
  
  
  return merged_clusters;
}







/* Another MergeClusteringResultsX function...
 *
 * Method (MergeClusteringResults2): Monti et al, http://www.cs.utexas.edu/users/ml/biodm/papers/MLJ-biodm2.pdf
 *
 */
ClusteringResult
MergeClusteringResults2( const vector<ClusteringResult> & crs, const int new_N_clusters )
{
  // Make sure there are multiple ClusteringResults and they all have the same number of contigs.
  assert( crs.size() >= 2 );
  int N_contigs = crs[0].NContigs();
  for ( size_t i = 1; i < crs.size(); i++ )
    assert( N_contigs == crs[i].NContigs() );
  
  // TODO: de-hardwire
  vector<double> weights;
  weights.push_back(1);
  weights.push_back(2);
  bool use_denominator = false;
  
  
  cout << Time() << ": Pre-processing" << endl;
  
  // Determine the number of unique "contig clustering codes".  Each code represents how a contig is clustered (or not) in each clustering,
  // The total number of possible codes is equal to product( N_i + 1 ), where N_i is the number of clusters in each clustering, and +1 is added because of the
  // possibility of a contig not being clustered in a particular clustering.
  vector<int> contig_clustering_code_offsets;
  int N_codes = 1;
  for ( size_t i = 0; i < crs.size(); i++ ) {
    contig_clustering_code_offsets.push_back( N_codes );
    PRINT3( i, crs[i]._N_clusters, N_codes );
    N_codes *= ( crs[i]._N_clusters + 1 );
    assert( N_codes > 0 ); // avoid integer overflow
  }
  
  // Find the clustering codes for each contig.
  vector<int> contig_clustering_code( N_contigs, 0 );
  for ( size_t i = 0; i < crs.size(); i++ )
    for ( int j = 0; j < N_contigs; j++ ) {
      int cID = crs[i]._contig_clusters[j];
      cID++; // increment by 1 so that -1 (not clustered) becomes 0 and cluster #X becomes X+1
      contig_clustering_code[j] += contig_clustering_code_offsets[i] * cID;
    }
  
  // As a pre-processing step (to save runtime later), make a mapping of clustering code to contig ID, and use it to condense contigs with the same clustering
  // code (i.e., they've always been placed into the same cluster) into a single node in the consensus matrix below.
  // Note that two contigs will not be condensed into a single node if the set of clusterings in which they are clustered is different.
  typedef multimap< int, int > MMAP;
  MMAP CCC_to_contig;
  
  for ( int i = 0; i < N_contigs; i++ )
    //PRINT2( i, contig_clustering_code[i] );
    CCC_to_contig.insert( make_pair( contig_clustering_code[i], i ) );
  
  
  
  vector<int> contig_to_node( N_contigs, -1 ); // vector of contigs mapping onto nodes created by having the same contig_clustering_code; -1 = not clustered
  vector<int> node_to_contig; // node representing a cluster of contigs -> a representative contig from that node
  vector<int> node_sizes;
  
  pair< MMAP::const_iterator, MMAP::const_iterator > iters;
  int key = CCC_to_contig.begin()->first;
  assert ( key == 0 );
  int N_counted = 0;
  
  // The following loop is like a regular iterator loop, except that in each iteration, it gets an entire equal_range of a multimap - i.e., the whole set of
  // key-value pairs for a given key.  The loop ends when iters.first = iters.second, which occurs when there are no more keys in the multimap.
  for ( iters = CCC_to_contig.equal_range( key );
	iters.first != iters.second;
	iters = CCC_to_contig.equal_range( key ) ) {
    
    // Create a node representing all of these contigs with the same contig clustering code.  Record one contig as a representative of this node.
    int node_ID = node_to_contig.size();
    int node_size = 0;
    for ( MMAP::const_iterator it = iters.first; it != iters.second; ++it ) {
      contig_to_node[ it->second ] = node_ID;
      N_counted++;
      node_size++;
    }
    
    node_to_contig.push_back( iters.first->second );
    node_sizes    .push_back( node_size );
    
    key = iters.second->first; // prepare for the next loop
  }
  
  
  PRINT( N_counted );
  
  
  
  // Make the "consensus matrix" of all clusterings, as defined in the Monti paper.
  // The rows and columns in the matrix are nodes (of contigs that are always clustered in the same way.)  Each element M(i,j) in the consensus matrix is
  // defined as the number of clusterings in which nodes i,j appear in the same cluster, divided by the total number of clusterings in which they are both
  // clustered.  (If the denominator is 0, then by definition, leave M(i,j) at 0.)
  int N_nodes = node_to_contig.size();
  
  
  cout << Time() << ": Allocating memory for consensus matrix with " << N_nodes << " rows/columns" << endl;
  boost::numeric::ublas::matrix<double> consensus_matrix( N_nodes, N_nodes, 0 );
  cout << Time() << ": Filling consensus matrix" << endl;
  
  
  
  // Loop over all pairs of nodes (i,j).  Note that we skip the case i,j=0, because node 0 has the contig clustering code 0, i.e., contigs never clustered.
  for ( int i = 1; i < N_nodes; i++ ) {
    for ( int j = i+1; j < N_nodes; j++ ) {
      
      double numerator = 0, denominator = 0;
      
      // Using the representative contigs to determine how these nodes are clustered, loop over each ClusteringResult and fill this element in the
      // consensus matrix.
      for ( size_t k = 0; k < crs.size(); k++ ) {
	
	int c1 = crs[k].cluster_ID( node_to_contig[i] );
	int c2 = crs[k].cluster_ID( node_to_contig[j] );
	
	if ( c1 == 0 || c2 == 0 ) continue; // ignore clusterings where the contigs are not both clustered
	
	denominator += weights[k];
	if ( c1 == c2 ) numerator += weights[k];
      }
      
      // If either the numerator or the denominator are zero, the matrix element is going to be zero.
      if ( numerator == 0 || denominator == 0 ) continue;
      
      // Make the matrix element.  Note that the consensus matrix, unlike the original link matrix, is symmetric instead of upper triangular.
      double Mij = numerator;
      if ( use_denominator ) Mij /= denominator;
      consensus_matrix(i,j) += Mij;
      consensus_matrix(j,i) += Mij;
      
    }
  }
  
  
  
  // Apply a simple form of hierarchical clustering to the consensus matrix.
  vector<int> node_clusters = AgglomerativeHierarchicalClusteringSimple( consensus_matrix, new_N_clusters );
  
  // These clustering results are by node.  Map them onto the contig IDs.
  vector<int> contig_clusters( N_contigs, -1 );
  for ( int i = 0 ; i < N_contigs; i++ )
    contig_clusters[i] = node_clusters[ contig_to_node[i] ];
  
  // Convert the hierarchical clustering result into a ClusteringResult object, and return it.
  return ClusteringResult( contig_clusters, crs[0]._contig_lengths );
}







/* Another MergeClusteringResultsX function...
 *
 * Method (MergeClusteringResults2): Dimitriadou et al, http://epub.wu.ac.at/94/1/document.pdf
 *
 */
ClusteringResult
MergeClusteringResults3( const vector<ClusteringResult> & crs, const int new_N_clusters )
{
  // Make sure there are multiple ClusteringResults and they all have the same number of contigs and the same number of clusters.
  assert( crs.size() >= 2 );
  int N_contigs = crs[0].NContigs();
  int N_clusters = crs[0].NClusters();
  for ( size_t i = 1; i < crs.size(); i++ ) {
    assert( N_contigs  == crs[i].NContigs() );
    assert( N_clusters == crs[1].NClusters() );
  }
  
  
  // Just work on the n=2 case for now.  This requires only one voting iteration.
  assert( crs.size() == 2 );
  
  
  vector<int> Ns_clusters; // number of clusters in each ClusteringResult
  for ( size_t i = 0; i < crs.size(); i++ )
    Ns_clusters.push_back( crs[i].NClusters() );
  
  
  // Apply the voting procedure.
  
  const ClusteringResult & C1 = crs[0];
  const ClusteringResult & C2 = crs[1];
  
  // For a given pair of clusterings C1, C2, find the length of all shared contigs between each cluster in C1 and each cluster in C2.
  cout << Time() << ": Finding shared contigs" << endl;
  vector< vector<int> > len_shared( N_clusters, vector<int>( N_clusters, 0 ) );
  
  // Loop over every contig and see how it's clustered in both clusterings.
  for ( int i = 0; i < N_contigs; i++ ) {
    int cluster1 = C1.cluster_ID(i);
    int cluster2 = C2.cluster_ID(i);
    if ( cluster1 == -1 ) continue;
    if ( cluster2 == -1 ) continue;
    
    //len_shared[cluster1][cluster2]++;
    len_shared[cluster1][cluster2] += C1._contig_lengths[i];
  }
  
  
  
  for ( int i = 0; i < N_clusters; i++ ) {
    cout << "len_shared[" << i << "]:\t";
    for ( int j = 0; j < N_clusters; j++ ) {
      cout << "\t" << len_shared[i][j];
    }
    cout << endl;
  }

  cout << Time() << ": Looping over shared contigs, determining cluster matches" << endl;
  
  vector<bool> used_C1( N_clusters, false );
  vector<bool> used_C2( N_clusters, false );
  
  // Find the best matching between the clusters in C1 and the clusters in C2.  This is done by a greedy algorithm.
  vector<int> match_C1C2( N_clusters, -1 );
  
  for ( int i = 0; i < N_clusters; i++ ) {
    
    int best_len_shared = -1;
    int best_j = -1, best_k = -1;
    
    // In each iteration, find the cluster in C1 and the cluster in C2 with the highest number of shared nodes.
    for ( int j = 0; j < N_clusters; j++ ) {
      if ( used_C1[j] ) continue;
      
      for ( int k = 0; k < N_clusters; k++ ) {
	if ( used_C2[k] ) continue;
	
	if ( len_shared[j][k] > best_len_shared ) {
	  best_len_shared = len_shared[j][k];
	  best_j = j;
	  best_k = k;
	}
      }
    }
    
    // Once this pair is found, record it, and then mark these two clusters to be ignored for future iterations.
    assert( best_j != -1 && best_k != -1 );
    
    match_C1C2[best_j] = best_k;
    PRINT4( i, best_j, best_k, best_len_shared );
    used_C1[best_j] = true;
    used_C2[best_k] = true;
    
  }
  
  
  
  // Create the matrix D2, which describes the merged result of C1 and C2.
  // Every contig clustered into C1 and C2 is assigned a cluster in D2 with a value of 1 (if the contig is clustered concordantly) or 0.5 (to both of its
  // clusters, if discordant.)
  // Note that D2 implicitly uses the cluster ordering in C2, not C1.  This will need to be fixed up if multi- (more than 2-) way merging is ever implemented.
  // Also find the "sureness" of each clustering.
  vector< vector<double> > D2( N_contigs, vector<double>( N_clusters, 0 ) );
  
  int N_clustered = 0, N_concordant = 0, N_discordant = 0;
  vector<double> sureness( N_contigs, 0 );
  
  for ( int i = 0; i < N_contigs; i++ ) {
    
    int cluster_in_C1 = C1.cluster_ID(i); // in C1 cluster ordering
    int cluster_in_C2 = C2.cluster_ID(i); // in C2 cluster ordering
    
    // If the contig wasn't clustered in either ordering, it's not clustered in D2, either.
    if ( cluster_in_C1 == -1 && cluster_in_C2 == -1 ) continue;
    
    N_clustered++;
    
    // If the contig was clustered in only one clustering, that one prevails.
    if ( cluster_in_C1 == -1 ) {
      //PRINT2( i, cluster_in_C2 );
      D2[i][ cluster_in_C2 ] = 1;
      sureness[i] = 1;
    }
    else if ( cluster_in_C2 == -1 ) {
      //PRINT2( i, match_C1C2[ cluster_in_C1 ] );
      D2[i][ match_C1C2[ cluster_in_C1 ] ] = 1;
      sureness[i] = 1;
    }
    
    // If the contig was clustered in both clusterings, each clustering assignment gets a score of 0.5.  Note that most contigs will be concordantly
    // clustered, so the 0.5's will fall in the same bin and add to 1.
    else {
      //PRINT3( i, match_C1C2[ cluster_in_C1 ], cluster_in_C2 );
      D2[i][ match_C1C2[ cluster_in_C1 ] ] += 0.5;
      D2[i][             cluster_in_C2   ] += 0.5;
      
      if ( match_C1C2[ cluster_in_C1 ] == cluster_in_C2 ) { sureness[i] = 1; N_concordant++; }
      else { sureness[i] = 0.5; N_discordant++; }
    }
    
  }
  
  double total_sureness = accumulate( sureness.begin(), sureness.end(), 0.0 );
  double avg_sureness = total_sureness / N_clustered;
  PRINT5( N_clustered, N_concordant, N_discordant, total_sureness, avg_sureness );
  
  
  
  // god damn it I just realized this is basically the same algorithm as in the Monti paper after wasting like an hour and a half of coding
  
  return crs[1];
}
