// For documentation, see HierarchicalClustering.h
#include "HierarchicalClustering.h"

// C libraries
#include <assert.h>

// STL declarations
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <iostream>
#include <algorithm> // count

#include "TimeMem.h"
#include "LinkMatrix.h" // LinkMatrixInt, LinkMatrixDouble
#include "TextFileParsers.h" // ParseTabDelimFile

#include "gtools/SAMStepper.h" // NTargetsInSAM, TargetCoverages

#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/matrix.hpp> // matrix












// FindShotgunCoverages: Return a vector indicating the coverage by shotgun reads of each contig in the draft genome.
// Note that the "coverage" here is defined as the NUMBER of reads aligning to each contig, per contig bp.  Read lengths are NOT taken into account.
// This function looks for a file at <shotgun_BAM>.covs.  If that file exists, it parses it, thus saving BAM-parsing time.  If not, it creates it.
vector<double>
FindShotgunCoverages( const string & shotgun_BAM, const int N_contigs )
{
  //PRINT( shotgun_BAM );
  
  // Firstly, if shotgun_BAM is an empty string, then the user has no information about shotgun coverage.
  // Return a vector of 1's, which will cause no link re-weighting.
  if ( shotgun_BAM == "" ) return vector<double>( N_contigs, 1 );
  
  // Make sure the BAM file matches this assembly.
  assert( boost::filesystem::is_regular_file( shotgun_BAM ) );
  assert( NTargetsInSAM( shotgun_BAM ) == N_contigs );
  
  // Look for the cache file at <shotgun_BAM>.covs.  If there's a file there, just read it and return the numbers in it.
  string cache_file = shotgun_BAM + ".covs";
  
  if ( boost::filesystem::is_regular_file( cache_file ) ) {
    vector<double> covs = ParseTabDelimFile<double>( cache_file, 0 );
    assert( (int) covs.size() == N_contigs );
    return covs;
  }
  
  
  // If there's no cache file, then we must parse the BAM file.  This is time-consuming.
  cout << Time() << ": Parsing shotgun BAM file " << shotgun_BAM << "..." << endl;
  vector<double> covs = TargetCoverages( shotgun_BAM );
  assert( (int) covs.size() == N_contigs );
  
  // Write the coverages to cache so we won't have to parse them from the BAM in the future.
  ofstream out( cache_file.c_str(), ios::out);
  for ( int i = 0; i < N_contigs; i++ )
    out << covs[i] << endl;
  out.close();
  
  return covs;
}






bool use_avg_link_density = false; // TEMP: tweak the linkage scores by the cluster densities



// TwoClusters: Helper struct for the AgglomerativeHierarchicalClustering function, below.
struct TwoClusters {
  int c1, c2;
  double link_weight; // normalized link weight between contigs in cluster 1 and contigs in cluster 2
  
  TwoClusters( const int & c1_in, const int & c2_in, const double & link_weight_in ) : c1(c1_in), c2(c2_in), link_weight(link_weight_in) {}
};




vector<int>
AgglomerativeHierarchicalClustering( const LinkMatrixInt & matrix_orig, const LinkMatrixDouble & matrix_in, const vector<int> & norms, const int N_clusters, const int MIN_CLUSTER_NORM )
{
  cout << Time() << ": HierarchicalClustering with N_clusters = " << N_clusters << endl;
  
  // Sanity checks.
  size_t N = matrix_in.size1();
  assert( N == matrix_in.size2() );
  assert( N == norms.size() );
  
  
  // Heuristic parameters that affect runtime/output only.
  static const int PRUNE_RATE = 10; // prune after each time we do 1/PRUNE_RATE of the total remaining number of merges
  // If two clusters both have this many contigs, keep a cache of their total linkages.
  // Lower numbers cause aggressive caching, which is faster - but only up to a point, and it's also more memory-intensive.
  static const int MIN_SIZE_FOR_CACHE = 20;
  static const bool DO_ENRICHMENT_CURVE = ( N_clusters == 1 );

  bool down_enrichment_curve = false;
  
  PRINT3( MIN_CLUSTER_NORM, MIN_SIZE_FOR_CACHE, DO_ENRICHMENT_CURVE );
  
  vector<bool> isolated_contigs = matrix_in._isolated_contigs;
  assert( N == isolated_contigs.size() );
  int N_non_isolated = count( isolated_contigs.begin(), isolated_contigs.end(), false );
  
  
  
  // These objects will contain intermediate products of the algorithm.
  vector<bool>   cluster_exists ( 2*N, false );
  vector<int>    cluster_size   ( 2*N, 0 ); // number of contigs in each cluster
  vector<int>    cluster_norm   ( 2*N, 0 ); // total length (norm) of the contigs of each cluster
  vector<double> cluster_N_links( 2*N, 0 ); // weight of Hi-C links among the contigs of each cluster
  vector<int>   contig_clusterID( N, -1 ); // final result
  
  // Put each contig in its own distinct cluster (except contigs marked with isolated_contig, which get left out for now).
  // This is the starting point of agglomerative clustering.
  for ( size_t i = 0; i < N; i++ ) {
    if ( isolated_contigs[i] ) continue;
    cluster_exists[i] = true;
    cluster_size[i] = 1;
    cluster_norm[i] = norms[i];
    cluster_N_links[i] = matrix_in(i,i); if ( cluster_N_links[i] == 0 ) cluster_N_links[i] = 1;
    contig_clusterID[i] = i;
  }
  
  
  // Calculate all possible "merge scores" for all pairs of clusters.  Make a list sorted by distance.
  // The initial "merge score" values for the initial (one-bin) clusters is simply the amount of link data between each pair of bins.
  // Also find the average link density between two contigs in the merge score map.
  cout << Time() << ": Creating a 'merge score map'..." << endl;
  typedef multimap< double, TwoClusters, greater<double> > MergeScoreMap;
  MergeScoreMap merge_score_map;
  double avg_link_density = 0;
  for ( LinkMatrixDouble::const_iterator1 it1 = matrix_in.begin1(); it1 != matrix_in.end1(); ++it1 )
    for ( LinkMatrixDouble::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2 )
      if ( !isolated_contigs[ it2.index1() ] && !isolated_contigs[ it2.index2() ] ) {
	merge_score_map.insert( make_pair( *it2, TwoClusters( it2.index1(), it2.index2(), *it2 ) ) );
	avg_link_density += *it2;
      }
  
  
  int N_possible_links = N_non_isolated * ( N_non_isolated - 1 ) / 2;
  assert( N_possible_links > 0 ); // if this fails, the input matrix is weird and empty (or you've set ISOLATED_COMPONENT_SIZE to 0)
  avg_link_density /= N_possible_links;
  if ( use_avg_link_density ) PRINT( avg_link_density );
  
  
  
  // Make a cache data structure to save runtime.
  // The data structure maps a pair of cluster numbers onto the number of links between the clusters (non-normalized).  Only add the info to this map for
  // clusters of size >50.  This is redundant with the information in merge_score_map, but it's much faster than recalculating it with each merge, especially
  // toward the end of clustering.
  map< pair<int,int>, double > linkage_cache;
  
  
  int N_merges = 0;
  int N_merges_since_prune = 0;
  int N_non_orig_clusters = 0;
  int N_large_clusters = 0; // "large" = total norm size >= MIN_CLUSTER_NORM
  bool past_start = false;
  
  
  // Hierarchical clustering!  Repeatedly perform the following three steps:
  // 1. Find the pair of clusters with the highest "merge score".
  // 2. Merge this pair and record the merging.
  // 3. If we've reached the optimal number of clusters, break.
  // 4. Update the multimap<> of merge scores to reflect the merging.
  cout << Time() << ": Doing the clustering!" << endl;
  PRINT2( N, N_non_isolated );
  
  
  while ( 1 ) {
    
    // 1. Find the pair of clusters with the highest "merge score".
    // This is easy - it's simply the first item in the map.
    MergeScoreMap::iterator it;
    
    for ( it = merge_score_map.begin(); it != merge_score_map.end(); ++it )
      if ( cluster_exists[ it->second.c1 ] && cluster_exists[ it->second.c2 ] )
	break;
    
    if ( it == merge_score_map.end() ) { cout << "empty merge score map. weird." << endl; break; }
    
    double best_linkage = it->first;
    assert( best_linkage > 0 );
    int best_i = it->second.c1;
    int best_j = it->second.c2;
    if ( best_i > best_j ) { int swap = best_i; best_i = best_j; best_j = swap; }
    
    
    // 2. Merge this pair and record the merging.
    
    // Merge together clusters best_i, best_j into a new cluster.
    // We do this instead of merging one cluster into another so that we don't have to remove all the existing entries in the merge score map that correspond
    // to the existing clusters; we can just mark these clusters as unused so that we can skip over those entries.
    assert( cluster_exists[best_i] );
    assert( cluster_exists[best_j] );
    
    // Get the new cluster ID.
    int new_cluster_ID = N + N_merges;
    assert( !cluster_exists[new_cluster_ID] );
    
    bool input_clusters_cached = ( cluster_size[best_i] >= MIN_SIZE_FOR_CACHE && cluster_size[best_j] >= MIN_SIZE_FOR_CACHE );
    
    // Record info for the new cluster.
    cluster_exists [new_cluster_ID] = true;
    cluster_size   [new_cluster_ID] = cluster_size   [best_i] + cluster_size   [best_j];
    cluster_norm   [new_cluster_ID] = cluster_norm   [best_i] + cluster_norm   [best_j];
    cluster_N_links[new_cluster_ID] = cluster_N_links[best_i] + cluster_N_links[best_j] + it->second.link_weight;
    
    // Update contig_clusterID, and create a new_cluster set.
    set<int> new_cluster;
    for ( size_t bin = 0; bin < N; bin++ )
      if ( contig_clusterID[bin] == best_i || contig_clusterID[bin] == best_j ) {
	contig_clusterID[bin] = new_cluster_ID;
	new_cluster.insert(bin);
      }
    
    
    // Keep tallies of the number of non-original (i.e., non-singleton) clusters and large clusters that have been created.
    N_non_orig_clusters++;
    if ( best_i >= int(N) ) N_non_orig_clusters--;
    if ( best_j >= int(N) ) N_non_orig_clusters--;
    if ( cluster_norm[new_cluster_ID] >= MIN_CLUSTER_NORM ) N_large_clusters++;
    if ( best_i >= int(N) && cluster_norm[best_i] >= MIN_CLUSTER_NORM ) N_large_clusters--;
    if ( best_j >= int(N) && cluster_norm[best_j] >= MIN_CLUSTER_NORM ) N_large_clusters--;
    
    
    //cout << "Merging clusters #" << best_i << " (" << cluster_norm[best_i] << "), #" << best_j << " (" << cluster_norm[best_j] << ") into #" << new_cluster_ID << endl;
    cout << Time() << ": After merge #" << N_merges << ": " << N_non_orig_clusters << " multi-node clusters, " << N_large_clusters << " large clusters, best linkage = " << best_linkage << endl;
    
    // If N_clusters is set to 1, then report on intra-cluster enrichment while clustering.  This is computationally expensive.
    if ( DO_ENRICHMENT_CURVE ) {
      if ( N_large_clusters >= 200 ) down_enrichment_curve = true;
      if ( down_enrichment_curve && N_large_clusters <= 100 )
      cout << "\tN LARGE CLUSTERS = " << N_large_clusters << "\tE(N) = " << matrix_orig.IntraClusterEnrichment( contig_clusterID, norms, cluster_norm, MIN_CLUSTER_NORM ) << endl;
    }
    
    N_merges++;
    N_merges_since_prune++;
    
    
    // Clear out info about the old clusters.
    cluster_exists[best_i] = false;
    cluster_exists[best_j] = false;
    cluster_size[best_i] = 0;
    cluster_size[best_j] = 0;
    cluster_norm[best_i] = 0;
    cluster_norm[best_j] = 0;
    cluster_N_links[best_i] = 0;
    cluster_N_links[best_j] = 0;
    
    
    
    
    // 3. If we've reached the optimal number of clusters, break.
    // Keep track of how many clusters are "non-original" - that is, not among the first N singleton clusters in the range [0,N).
    // This number will start at 0, then increase as clustering progresses, then decrease.  We want to catch it as it hits N_clusters on the way down.
    // This is better than counting the total number of clusters, because a lot of contigs have no data and will always remain in their original clusters.
    
    //if ( N_non_orig_clusters > N_clusters * 2 ) past_start = true; // the *2 provides a buffer against small fluctuations
    if ( N_large_clusters > N_clusters * 2 ) past_start = true; // the *2 provides a buffer against small fluctuations
    
    
    if ( past_start && N_large_clusters == N_clusters ) {
      cout << N_merges << " merges made so far; this leaves " << N_clusters << " clusters, and so we're done!" << endl;
      break;
    }
    
    
    
    
    
    // 4. Update the multimap<> of merge scores to reflect the merging.  This has two sub-parts, 4a and 4b.
    
    
    // 4a. Find all merge scores involving clusters that no longer exist, and remove them.
    // We only do this rarely because it's very time-consuming.
    if ( N_merges_since_prune > int ( N_non_isolated - N_merges ) / PRUNE_RATE ) {
      N_merges_since_prune = 0;
      
      for ( it = merge_score_map.begin(); it != merge_score_map.end(); it++ ) {
	while ( it != merge_score_map.end() &&
		( !cluster_exists[ it->second.c1 ] ||
		  !cluster_exists[ it->second.c2 ] ) )
	  merge_score_map.erase(it++);
	
	if ( it == merge_score_map.end() ) break;
      }
    }
    
    
    
    // 4b. Calculate new score entries for the new cluster - that is, the average linkage from this cluster to each other cluster.
    // TODO: this runtime fix works! Propagate it to Lachesis
    vector<double> total_linkage_by_cluster( 2*N, 0 );
    
    for ( int i = 0; i < int(N); i++ ) {
      int cluster_ID = contig_clusterID[i];
      //PRINT2( i, contig_clusterID[i] );
      if ( cluster_ID == new_cluster_ID ) continue; // no need to calculate linkages within a cluster
      if ( cluster_ID == -1 ) continue; // this happens if contig_isolated[i]
      
      // If both of the constituent clusters of new_cluster are sufficiently large, and the loop cluster is also sufficiently large, the linkages should
      // already be cached.
      if ( input_clusters_cached && cluster_size[cluster_ID] >= MIN_SIZE_FOR_CACHE ) {
	if ( total_linkage_by_cluster[cluster_ID] != 0 ) continue; // already saw a contig from this cluster, no need to recalculate
	map< pair<int,int>, double >::const_iterator it1, it2;
	it1 = linkage_cache.find( make_pair( min( cluster_ID, best_i ), max( cluster_ID, best_i ) ) );
	it2 = linkage_cache.find( make_pair( min( cluster_ID, best_j ), max( cluster_ID, best_j ) ) );
	assert( it1 != linkage_cache.end() );
	assert( it2 != linkage_cache.end() );
	double linkage = it1->second + it2->second;
	//PRINT6( cluster_ID, new_cluster_ID, cluster_size[cluster_ID], cluster_size[new_cluster_ID], linkage );
	linkage_cache.insert( make_pair( make_pair( min( cluster_ID, new_cluster_ID ), max( cluster_ID, new_cluster_ID ) ), linkage ) );
	total_linkage_by_cluster[cluster_ID] = linkage;
      }
      
      else // no cached data (normal case)
	for ( set<int>::const_iterator it = new_cluster.begin(); it != new_cluster.end(); ++it )
	  total_linkage_by_cluster[cluster_ID] += matrix_in( min(i,*it),max(i,*it) );
    }
    
    
    
    
    // If this cluster has a lot of nodes, find its linkages with other clusters and add this info to the linkage cache.
    if ( cluster_size[new_cluster_ID] >= MIN_SIZE_FOR_CACHE )
      for ( int i = 0; i < int(2*N); i++ ) {
	if ( i == new_cluster_ID ) continue;
	if ( cluster_size[i] >= MIN_SIZE_FOR_CACHE ) {
	  //cout << "Adding linkage_cache[" << min( i, new_cluster_ID ) << "," << max( i, new_cluster_ID ) << "]" << endl;
	  linkage_cache.insert( make_pair( make_pair( min( i, new_cluster_ID ), max( i, new_cluster_ID ) ), total_linkage_by_cluster[i] ) );
	}
      }
    
    
    
    
    
    for ( int i = 0; i < int(2*N); i++ ) {
      if ( total_linkage_by_cluster[i] == 0 ) continue;
      
      assert( cluster_exists[i] );
      double avg_linkage = double( total_linkage_by_cluster[i] ) / cluster_size[i] / cluster_size[new_cluster_ID];
      if ( use_avg_link_density ) {
	//cout << "avg linkage before = " << avg_linkage << endl;
	//PRINT2( cluster_N_links[i], cluster_N_links[new_cluster_ID] );
	double cluster_expected_N_links_i = cluster_norm[i]              * cluster_norm[i]              * avg_link_density;
	double cluster_expected_N_links_n = cluster_norm[new_cluster_ID] * cluster_norm[new_cluster_ID] * avg_link_density;
	double density_norm = ( cluster_N_links[i] / cluster_expected_N_links_i ) * ( cluster_N_links[new_cluster_ID] / cluster_expected_N_links_n );
	density_norm *= ( cluster_size[i] * cluster_size[new_cluster_ID] );
	PRINT4( i, new_cluster_ID, avg_linkage, density_norm );
	avg_linkage /= density_norm;
	//cout << "avg linkage after  = " << avg_linkage << endl;
      }
      //cout << "Adding a link between:\t" << i << "\t" << new_cluster_ID << " with linkage " << avg_linkage << endl;
      assert( i < new_cluster_ID );
      merge_score_map.insert( make_pair( avg_linkage, TwoClusters( i, new_cluster_ID, total_linkage_by_cluster[i] ) ) );
    }
    
    
  }
  
  
  cout << Time() << ": Done clustering!" << endl;
  
  // The cluster assignments are now in the variable contig_clusterID, as follows:
  // CASE 1: contig_clusterID[i] = cID, cID >= N: Contig has been placed into cluster #cID.  Cluster #cID might be too small, though (Case 1a).
  // CASE 2: contig_clusterID[i] = -1: Contig was marked as isolated, thus ineligible for clustering.
  // CASE 3: contig_clusterID[i] = i: Contig was eligible for clustering but was left as a singleton due to incomplete clustering.
  // Find cases 1a and 3, and make them -1!
  int N_eligible_contigs_unclustered = 0;
  for ( int i = 0; i < int(N); i++ )
    if ( contig_clusterID[i] >= 0 && contig_clusterID[i] < int(N) ) { // case 1a
      assert( contig_clusterID[i] == i );
      cluster_size[i] = 0;
      contig_clusterID[i] = -1;
      N_eligible_contigs_unclustered++;
    }
    else if ( cluster_norm[ contig_clusterID[i] ] < MIN_CLUSTER_NORM ) { // case 3
      cluster_size[i] = 0;
      contig_clusterID[i] = -1;
      N_eligible_contigs_unclustered++;
    }
  
  //PRINT( N_eligible_contigs_unclustered );
  
  //for ( size_t i = N; i < 2*N; i++ )
  //if ( cluster_size[i] > 1 )
  //PRINT4( i, cluster_size[i], cluster_norm[i], cluster_N_links[i] );
  
  
  
  // Sort the clusters!
  // Count the size of each cluster, and use these numbers to renumber the clusters so that cluster #0 is the largest, etc.
  multimap< int, int, greater<int> > clusters_sorted_by_size;
  for ( size_t i = N; i < 2*N; i++ )
    clusters_sorted_by_size.insert( make_pair( cluster_size[i], i ) );
  
  map<int,int> cluster_renumbers;
  cluster_renumbers[-1] = -1; // un-clustered contigs should stay un-clustered under renumbering
  int renumber_ID = 0;
  for ( multimap< int, int, greater<int> >::const_iterator it = clusters_sorted_by_size.begin(); it != clusters_sorted_by_size.end(); ++it ) {
    //cout << "Cluster #" << it->second << " (renumbered to " << renumber_ID << ") has size " << it->first << "\n";
    cluster_renumbers[ it->second ] = renumber_ID++;
  }
  
  
  vector<int> contig_clusterID_renumbered( N, 0 );
  for ( size_t i = 0; i < N; i++ )
    contig_clusterID_renumbered[i] = cluster_renumbers[ contig_clusterID[i] ];
  
  
  return contig_clusterID_renumbered;
}













vector<int>
AgglomerativeHierarchicalClusteringSimple( const boost::numeric::ublas::matrix<double> & matrix_in, const int N_clusters )
{
  cout << Time() << ": HierarchicalClusteringSimple with N_clusters = " << N_clusters << endl;
  
  // Sanity checks.
  size_t N = matrix_in.size1();
  assert( N == matrix_in.size2() );
  
  
  // Heuristic parameters that affect runtime only.
  static const int PRUNE_RATE = 10; // prune after each time we do 1/PRUNE_RATE of the total remaining number of merges
  // If two clusters both have this many contigs, keep a cache of their total linkages.
  // Lower numbers cause aggressive caching, which is faster - but only up to a point, and it's also more memory-intensive.
  static const int MIN_SIZE_FOR_CACHE = 20;
  
  
  
  
  // These objects will contain intermediate products of the algorithm.
  vector<bool>   cluster_exists ( 2*N, false );
  vector<int>    cluster_size   ( 2*N, 0 ); // number of contigs in each cluster
  vector<int>   contig_clusterID( N, -1 ); // final result
  
  // Put each contig in its own distinct cluster (except contigs marked with isolated_contig, which get left out for now).
  // This is the starting point of agglomerative clustering.
  for ( size_t i = 0; i < N; i++ ) {
    cluster_exists[i] = true;
    cluster_size[i] = 1;
    contig_clusterID[i] = i;
  }
  
  
  // Calculate all possible "merge scores" for all pairs of clusters.  Make a list sorted by distance.
  // The initial "merge score" values for the initial (one-bin) clusters is simply the amount of link data between each pair of bins.
  // Also find the average link density between two contigs in the merge score map.
  cout << Time() << ": Creating a 'merge score map'..." << endl;
  typedef multimap< double, TwoClusters, greater<double> > MergeScoreMap;
  MergeScoreMap merge_score_map;
  for ( boost::numeric::ublas::matrix<double>::const_iterator1 it1 = matrix_in.begin1(); it1 != matrix_in.end1(); ++it1 )
    for ( boost::numeric::ublas::matrix<double>::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2 )
      if ( it2.index1() != it2.index2() ) // skip on-diagonal elements
	merge_score_map.insert( make_pair( *it2, TwoClusters( it2.index1(), it2.index2(), *it2 ) ) );
  
  
  
  
  
  // Make a cache data structure to save runtime.
  // The data structure maps a pair of cluster numbers onto the number of links between the clusters (non-normalized).  Only add the info to this map for
  // clusters of size >50.  This is redundant with the information in merge_score_map, but it's much faster than recalculating it with each merge, especially
  // toward the end of clustering.
  map< pair<int,int>, double > linkage_cache;
  
  
  int N_merges = 0;
  int N_merges_since_prune = 0;
  int N_non_orig_clusters = 0;
  bool past_start = false;
  
  
  // Hierarchical clustering!  Repeatedly perform the following three steps:
  // 1. Find the pair of clusters with the highest "merge score".
  // 2. Merge this pair and record the merging.
  // 3. If we've reached the optimal number of clusters, break.
  // 4. Update the multimap<> of merge scores to reflect the merging.
  cout << Time() << ": Doing the clustering!" << endl;
  
  
  while ( 1 ) {
    
    // 1. Find the pair of clusters with the highest "merge score".
    // This is easy - it's simply the first item in the map.
    MergeScoreMap::iterator it;
    
    for ( it = merge_score_map.begin(); it != merge_score_map.end(); ++it )
      if ( cluster_exists[ it->second.c1 ] && cluster_exists[ it->second.c2 ] )
	break;
    
    if ( it == merge_score_map.end() ) { cout << "empty merge score map. weird." << endl; break; }
    
    double best_linkage = it->first;
    if ( best_linkage <= 0 ) { cout << "best_linkage = " << best_linkage << " <= 0.  weird." << endl; break; }
    assert( best_linkage > 0 );
    int best_i = it->second.c1;
    int best_j = it->second.c2;
    if ( best_i > best_j ) { int swap = best_i; best_i = best_j; best_j = swap; }
    
    
    // 2. Merge this pair and record the merging.
    
    // Merge together clusters best_i, best_j into a new cluster.
    // We do this instead of merging one cluster into another so that we don't have to remove all the existing entries in the merge score map that correspond
    // to the existing clusters; we can just mark these clusters as unused so that we can skip over those entries.
    assert( cluster_exists[best_i] );
    assert( cluster_exists[best_j] );
    
    // Get the new cluster ID.
    int new_cluster_ID = N + N_merges;
    assert( !cluster_exists[new_cluster_ID] );
    
    bool input_clusters_cached = ( cluster_size[best_i] >= MIN_SIZE_FOR_CACHE && cluster_size[best_j] >= MIN_SIZE_FOR_CACHE );
    
    // Record info for the new cluster.
    cluster_exists [new_cluster_ID] = true;
    cluster_size   [new_cluster_ID] = cluster_size   [best_i] + cluster_size   [best_j];
    
    // Update contig_clusterID, and create a new_cluster set.
    set<int> new_cluster;
    for ( size_t bin = 0; bin < N; bin++ )
      if ( contig_clusterID[bin] == best_i || contig_clusterID[bin] == best_j ) {
	contig_clusterID[bin] = new_cluster_ID;
	new_cluster.insert(bin);
      }
    
    
    // Keep tallies of the number of non-original (i.e., non-singleton) clusters and large clusters that have been created.
    N_non_orig_clusters++;
    if ( best_i >= int(N) ) N_non_orig_clusters--;
    if ( best_j >= int(N) ) N_non_orig_clusters--;
    
    
    //cout << "Merging clusters #" << best_i << ", #" << best_j << " into #" << new_cluster_ID << endl;
    cout << Time() << ": After merge #" << N_merges << ": " << N_non_orig_clusters << " multi-node clusters, best linkage = " << best_linkage << endl;
    
    N_merges++;
    N_merges_since_prune++;
    
    
    // Clear out info about the old clusters.
    cluster_exists[best_i] = false;
    cluster_exists[best_j] = false;
    cluster_size[best_i] = 0;
    cluster_size[best_j] = 0;
    
    
    
    
    // 3. If we've reached the optimal number of clusters, break.
    // Keep track of how many clusters are "non-original" - that is, not among the first N singleton clusters in the range [0,N).
    // This number will start at 0, then increase as clustering progresses, then decrease.  We want to catch it as it hits N_clusters on the way down.
    // This is better than counting the total number of clusters, because a lot of contigs have no data and will always remain in their original clusters.
    
    if ( N_non_orig_clusters > N_clusters + 1 ) past_start = true; // the +1 provides a buffer against small fluctuations
    
    //PRINT3( past_start, N_non_orig_clusters, N_clusters );
    if ( past_start && N_non_orig_clusters == N_clusters ) {
      cout << N_merges << " merges made so far; this leaves " << N_clusters << " clusters, and so we're done!" << endl;
      break;
    }
    
    
    
    
    
    // 4. Update the multimap<> of merge scores to reflect the merging.  This has two sub-parts, 4a and 4b.
    
    
    // 4a. Find all merge scores involving clusters that no longer exist, and remove them.
    // We only do this rarely because it's very time-consuming.
    if ( N_merges_since_prune > int ( N - N_merges ) / PRUNE_RATE ) {
      N_merges_since_prune = 0;
      
      for ( it = merge_score_map.begin(); it != merge_score_map.end(); it++ ) {
	while ( it != merge_score_map.end() &&
		( !cluster_exists[ it->second.c1 ] ||
		  !cluster_exists[ it->second.c2 ] ) )
	  merge_score_map.erase(it++);
	
	if ( it == merge_score_map.end() ) break;
      }
    }
    
    
    
    // 4b. Calculate new score entries for the new cluster - that is, the average linkage from this cluster to each other cluster.
    // TODO: this runtime fix works! Propagate it to Lachesis
    vector<double> total_linkage_by_cluster( 2*N, 0 );
    
    for ( int i = 0; i < int(N); i++ ) {
      int cluster_ID = contig_clusterID[i];
      //PRINT2( i, contig_clusterID[i] );
      if ( cluster_ID == new_cluster_ID ) continue; // no need to calculate linkages within a cluster
      if ( cluster_ID == -1 ) continue; // this happens if contig_isolated[i]
      
      // If both of the constituent clusters of new_cluster are sufficiently large, and the loop cluster is also sufficiently large, the linkages should
      // already be cached.
      if ( input_clusters_cached && cluster_size[cluster_ID] >= MIN_SIZE_FOR_CACHE ) {
	if ( total_linkage_by_cluster[cluster_ID] != 0 ) continue; // already saw a contig from this cluster, no need to recalculate
	map< pair<int,int>, double >::const_iterator it1, it2;
	it1 = linkage_cache.find( make_pair( min( cluster_ID, best_i ), max( cluster_ID, best_i ) ) );
	it2 = linkage_cache.find( make_pair( min( cluster_ID, best_j ), max( cluster_ID, best_j ) ) );
	assert( it1 != linkage_cache.end() );
	assert( it2 != linkage_cache.end() );
	double linkage = it1->second + it2->second;
	//PRINT6( cluster_ID, new_cluster_ID, cluster_size[cluster_ID], cluster_size[new_cluster_ID], linkage );
	linkage_cache.insert( make_pair( make_pair( min( cluster_ID, new_cluster_ID ), max( cluster_ID, new_cluster_ID ) ), linkage ) );
	total_linkage_by_cluster[cluster_ID] = linkage;
      }
      
      else // no cached data (normal case)
	for ( set<int>::const_iterator it = new_cluster.begin(); it != new_cluster.end(); ++it )
	  total_linkage_by_cluster[cluster_ID] += matrix_in( min(i,*it),max(i,*it) );
    }
    
    
    
    
    // If this cluster has a lot of nodes, find its linkages with other clusters and add this info to the linkage cache.
    if ( cluster_size[new_cluster_ID] >= MIN_SIZE_FOR_CACHE )
      for ( int i = 0; i < int(2*N); i++ ) {
	if ( i == new_cluster_ID ) continue;
	if ( cluster_size[i] >= MIN_SIZE_FOR_CACHE ) {
	  //cout << "Adding linkage_cache[" << min( i, new_cluster_ID ) << "," << max( i, new_cluster_ID ) << "]" << endl;
	  linkage_cache.insert( make_pair( make_pair( min( i, new_cluster_ID ), max( i, new_cluster_ID ) ), total_linkage_by_cluster[i] ) );
	}
      }
    
    
    
    
    
    for ( int i = 0; i < int(2*N); i++ ) {
      if ( total_linkage_by_cluster[i] == 0 ) continue;
      
      assert( cluster_exists[i] );
      double avg_linkage = double( total_linkage_by_cluster[i] ) / cluster_size[i] / cluster_size[new_cluster_ID];
      //cout << "Adding a link between:\t" << i << "\t" << new_cluster_ID << " with linkage " << avg_linkage << endl;
      assert( i < new_cluster_ID );
      merge_score_map.insert( make_pair( avg_linkage, TwoClusters( i, new_cluster_ID, total_linkage_by_cluster[i] ) ) );
    }
    
    
  }
  
  
  cout << Time() << ": Done clustering!" << endl;
  
  // The cluster assignments are now in the variable contig_clusterID, as follows:
  // contig_clusterID[i] = cID, cID >= N: Contig has been placed into cluster #cID.
  // contig_clusterID[i] = -1: Contig was marked as isolated, thus ineligible for clustering.
  // contig_clusterID[i] = i: Contig was eligible for clustering but was left as a singleton due to incomplete clustering.  Find these cases and make them -1!
  int N_eligible_contigs_unclustered = 0;
  for ( int i = 0; i < int(N); i++ )
    if ( contig_clusterID[i] >= 0 && contig_clusterID[i] < int(N) ) {
      assert( contig_clusterID[i] == i );
      contig_clusterID[i] = -1;
      N_eligible_contigs_unclustered++;
    }
  
  //PRINT( N_eligible_contigs_unclustered );
  
  
  
  
  // Sort the clusters!  (TODO: move this logic into the ClusteringResult constructor; allow multiple ways of sorting clusters)
  // Count the size of each (non-original) cluster, and use these numbers to renumber the clusters so that cluster #0 is the largest, etc.
  multimap< int, int, greater<int> > clusters_sorted_by_size;
  for ( size_t i = N; i < 2*N; i++ )
    if ( cluster_size[i] > 1 )
      clusters_sorted_by_size.insert( make_pair( cluster_size[i], i ) );
  
  map<int,int> cluster_renumbers;
  cluster_renumbers[-1] = -1; // un-clustered contigs should stay un-clustered under renumbering
  int renumber_ID = 0;
  for ( multimap< int, int, greater<int> >::const_iterator it = clusters_sorted_by_size.begin(); it != clusters_sorted_by_size.end(); ++it ) {
    //cout << "Cluster #" << it->second << " (renumbered to " << renumber_ID << ") has size " << it->first << "\n";
    cluster_renumbers[ it->second ] = renumber_ID++;
  }
  
  
  vector<int> contig_clusterID_renumbered( N, 0 );
  for ( size_t i = 0; i < N; i++ )
    contig_clusterID_renumbered[i] = cluster_renumbers[ contig_clusterID[i] ];
  
  
  return contig_clusterID_renumbered;
}

