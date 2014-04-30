// For documentation, see LinkMatrix.h
#include "LinkMatrix.h"

#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm> // sort, count
#include <stdlib.h> // srand48, drand48


#include "TimeMem.h"
#include "TrueMapping.h"
#include "TextFileParsers.h"

// Boost libraries
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>


// Template instantiations for the functions below.
template void   LinkMatrix<int64_t>::Resize( size_t );
template void   LinkMatrix<double >::Resize( size_t );
template void   LinkMatrix<int64_t>::Bootstrap( const bool );
template void   LinkMatrix<int64_t>::Normalize( const vector<int> &, LinkMatrix<double> & ) const;
template void   LinkMatrix<int64_t>::LimitNeighbors( const int, const bool );
template void   LinkMatrix<double >::LimitNeighbors( const int, const bool );
template void   LinkMatrix<int64_t>::FindIsolatedNodes( const int );
template void   LinkMatrix<double >::FindIsolatedNodes( const int );
template void   LinkMatrix<int64_t>::IsolateNormalizeLimit( const vector<int> &, const int, const int, const string &, LinkMatrix<double> & ) const;
template double LinkMatrix<int64_t>::IntraClusterEnrichment( const vector<int> &, const vector<int> &, const vector<int> &, const int ) const;
template void   LinkMatrix<int64_t>::PrintSparse() const;
template void   LinkMatrix<double> ::WriteXGMML( const string &, const double &, const TrueMapping * ) const;


template<class T> void LinkMatrix<T>::Resize( size_t N )
{
  this->resize( N, N, false );
  N_links = 0;
  N_pairs = 0;
  _isolated_contigs.clear();
}







// Bootstrap the Hi-C link data.  That is, re-create the dataset of links by using random sampling _with replacement_ from the current dataset.
// For proper bootstrap sampling, this function should be run once, followed by clustering, and the clustering results should be compared for >1k bootstrapping
// instances.
template<class T> void
LinkMatrix<T>::Bootstrap( const bool set_random_seed )
{
  cout << Time() << ": Bootstrapping the data..." << endl;
  
  // Make a mapping of random numbers in the range [0,N_pairs) to data points.  This will be used in determining the resample.
  // This is the most time- and memory-consuming step.
  vector< pair<int,int> > resample( N_pairs );
  int resample_N_pairs = 0;
  
  for ( typename LinkMatrix<T>::const_iterator1 it1 = this->begin1(); it1 != this->end1(); ++it1 )
    for ( typename LinkMatrix<T>::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2 ) {
      fill_n( resample.begin() + resample_N_pairs, *it2, make_pair( it2.index1(), it2.index2() ) );
      resample_N_pairs += *it2;
    }
  
  assert( resample_N_pairs == N_pairs ); // if this fails, the dataset size isn't actually _N_pairs
  
  
  // Clear the matrix in preparation for resampling.
  int N_links_orig = N_links;
  int N_pairs_orig = N_pairs;
  vector<bool> isolated_contigs = _isolated_contigs;
  Resize( size() );
  _isolated_contigs = isolated_contigs;
  N_pairs = N_pairs_orig;
  
  
  // Set the random seed, if desired.
  if ( set_random_seed ) srand48( time(0) );
  
  
  // Do the resampling!
  for ( int i = 0; i < N_pairs; i++ ) {
    int draw = lrand48() % N_pairs;
    const int & ID1 = resample[draw].first, & ID2 = resample[draw].second;
    if ( (*this)(ID1,ID2) == 0 ) N_links++;
    (*this)(ID1,ID2) += 1;
  }
  
  
  cout << Time() << ": Done bootstrapping!  The total number of pairs is constant at " << N_pairs << ", but the number of links (nonzero matrix entries) has decreased from " << N_links_orig << " to " << N_links << '.' << endl;
}



// Normalize: Normalize the link density in a matrix by contig norms (lengths).  Note that output is a LinkMatrix<double>.
template<typename int64_t> void
LinkMatrix<int64_t>::Normalize( const vector<int> & norms, LinkMatrix<double> & normed_matrix ) const
{
  size_t N = norms.size();
  assert( N == size() );
  
  // Find the magnitude of the average norm.  This is used to re-scale the numbers after normalization to make sure they stay in a tractable range.
  int magnitude_norm;
  {
    int magnitude_sum = 0;
    int magnitude_N = 0;
    for ( size_t i = 0; i < N; i++ ) {
      magnitude_sum += int( log10( norms[i] ) );
      magnitude_N++;
    }
    magnitude_norm = pow( 10, ( magnitude_sum / magnitude_N ) );
    //PRINT( magnitude_norm );
  }
  
  
  // Normalize the matrix data.
  normed_matrix.Resize( N );
  normed_matrix._isolated_contigs = _isolated_contigs;
  double sum = 0;
  
  for ( LinkMatrixInt::const_iterator1 it1 = this->begin1(); it1 != this->end1(); ++it1 )
    for ( LinkMatrixInt::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2 ) {
      if ( *it2 == 0 ) continue;
      int i = it2.index1(), j = it2.index2();
      double value = double( *it2 ) * magnitude_norm / norms[i] * magnitude_norm / norms[j];
      normed_matrix(i,j) = value;
      sum += value;
    }
  
  normed_matrix.N_links = N_links;
  normed_matrix.N_pairs = sum;
}





// LimitNeighbors: Remove links from the matrix so that every row has a maximum of K nonzero elements in it.  Specifically, remove all links L(a,b) from
// the graph unless L is one of the top K highest numbers in both its row and its column: that is, it's one of the top K hits for both nodes a and b.
// If reweight = true, also set the weights of all (non-removed) links to indicate the strength of overlapping nearest neighbors.
// This step is also called "K-nearest neighbor sparsification" and is an important pre-processing step in the Jarvis-Patrick and SNN clustering algorithms.
// See http://www-users.cs.umn.edu/~kumar/papers/kdd02_snn_28.pdf
template<class T> void
LinkMatrix<T>::LimitNeighbors( const int K, const bool reweight )
{
  assert( K > 0 );
  assert( K < size() );
  
  typedef multimap< T, int, greater<T> > ValueMap;
  typename ValueMap::const_iterator it1, it2;
  
  
  // For each node, compile a list of the nodes it's connected to and the strengths of those connections.  The numbers go in a ValueMap to sort them by value.
  PRINT( MemUsage() );
  vector<ValueMap> value_maps( size() );
  vector< vector<int> > NN_list( size(), vector<int>(K, 0) );
  
  for ( typename LinkMatrix<T>::const_iterator1 it1 = this->begin1(); it1 != this->end1(); ++it1 )
    for ( typename LinkMatrix<T>::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2 )
      if ( *it2 != 0 ) {
	value_maps[ it2.index1() ].insert( make_pair( *it2, it2.index2() ) );
	value_maps[ it2.index2() ].insert( make_pair( *it2, it2.index1() ) );
      }
  PRINT( MemUsage() );
  
  // Find all pairs of nodes i, j such that i is one of the K nearest neighbors of j, and vice versa.
  vector< pair<int,int> > good_pairs;
  for ( int i = 0; i < size(); i++ ) {
    
    // Consider the K nodes linked most strongly to i.
    int k1 = 0;
    for ( it1 = value_maps[i].begin(); it1 != value_maps[i].end() && k1 < K; ++it1, ++k1 ) {
      
      const int & j = it1->second;
      NN_list[i][k1] = j;
      if ( i > j ) continue; // if i > j, then we'll find this pair later when they're switched
      
      // Consider the K nodes linked most strongly to j.  If i is in this list, then we've found a pair.
      int k2 = 0;
      for ( it2 = value_maps[j].begin(); it2 != value_maps[j].end() && k2 < K; ++it2, ++k2 ) {
	if ( it2->second == i ) { // found a pair!
	  good_pairs.push_back( make_pair( i, j ) );
	  break;
	}
      }
      
    }
  }
  
  
  
  
  // Remake the matrix using only the connections that are among the top K connections for each of the nodes they connect.
  // Also, if reweight = true, reweight the connections in accordance with the number of shared neighbors the nodes have.
  cout << Time() << ": Reducing number of matrix elements from " << N_links << " to " << good_pairs.size() << endl;
  
  LinkMatrix<T> new_matrix;
  new_matrix.Resize( size() );
  new_matrix._isolated_contigs = _isolated_contigs;
  
  //for ( set< pair<int,int> >::const_iterator it = good2.begin(); it != good2.end(); ++it ) {
  //int i = it->first, j = it->second;
  for ( int p = 0; p < good_pairs.size(); p++ ) {
    int i = good_pairs[p].first, j = good_pairs[p].second;
    
    if ( !reweight ) { new_matrix(i,j) = (*this)(i,j); }
    else {
      
      // Apply the formula from the paper: strength( i, j ) = sum over m,n in [0,K-1] of (k-m) * (k-n), where i_m = j_m
      int new_value = K; // TODO: K != 0; weights should always be nonzero or some of them will be dropped; edges could have each other in their NN lists but not any others in common - this is a departure from the strict formula of the paper
      
      for ( int k1 = 0; k1 < K; k1++ )
	for ( int k2 = 0; k2 < K; k2++ )
	  if ( NN_list[i][k1] == NN_list[j][k2] )
	    new_value += (K-k1) * (K-k2);
      
      new_matrix(i,j) = new_value;
    }
    
    new_matrix.N_links++;
    new_matrix.N_pairs += new_matrix(i,j);
  }
  
  *this = new_matrix;
}





// FindIsolatedNodes: Fill _isolated_contigs, which indicating for each contig whether it's part of a component of size <= max_isolated_component_size.
template<class T> void
LinkMatrix<T>::FindIsolatedNodes( const int max_isolated_component_size )
{
  int N = size();
  
  // Convert the matrix into a Boost graph.
  typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS, boost::no_property > Graph;
  Graph G( N );
  for ( typename LinkMatrix<T>::const_iterator1 it1 = this->begin1(); it1 != this->end1(); ++it1 )
    for ( typename LinkMatrix<T>::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2 )
      add_edge( it2.index1(), it2.index2(), G );
  
  
  // Find the connected components in this graph.
  vector<int> CCs( N );
  int N_components = connected_components( G, &CCs[0] );
  
  
  // Find the number of contigs in each connected component.
  vector<int> N_contigs_per_component( N_components, 0 );
  for ( size_t i = 0; i < N; i++ )
    N_contigs_per_component[ CCs[i] ]++;
  
  
  // Determine each contig to be isolated if it's in a component that's too small.
  _isolated_contigs.resize( N, false );
  for ( size_t i = 0; i < N; i++ )
    if ( N_contigs_per_component[ CCs[i] ] <= max_isolated_component_size )
      _isolated_contigs[i] = true;
  
}



// IsolateNormalizeLimit: Use the results of FindIsolatedNodes(), Normalize(), and LimitNeighbors() in that order.
// Disconnect all nodes that are part of a component of size <= max_isolated_component_size; then normalize by contig norms; then run LimitNeighbors with K.
// shotgun_BAM: A BAM file of shotgun reads aligned to the draft assembly.  Used to adjust link weights, on the assumption that contigs with very different
//              coverages are less likely to belong to the same species cluster.  Can be an empty string, in which case no adjustment is done.
template<class T> void
LinkMatrix<T>::IsolateNormalizeLimit( const vector<int> & norms, const int max_isolated_component_size, const int K, const string & shotgun_BAM, LinkMatrix<double> & normed_matrix ) const
{
  PRINT2( max_isolated_component_size, K );
  
  /*
  // TEMP: Remove all single-pair links from the LinkMatrix.
  int MIN_N_PAIRS = 0;
  cout << Time() << ": Removing links of weight < " << MIN_N_PAIRS << endl;
  LinkMatrixInt matrix2;
  matrix2.Resize(N);
  matrix2.isolated_contigs = _isolated_contigs;
  for ( LinkMatrixInt::const_iterator1 it1 = matrix_in.begin1(); it1 != matrix_in.end1(); ++it1 )
    for ( LinkMatrixInt::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2 ) {
      if ( *it2 >= MIN_N_PAIRS )
	matrix2( it2.index1(), it2.index2() ) = *it2;
    }
  */
  
  
  
  
  // Normalize the matrix elements by contig lengths.  The result of this opetation goes into the normed_matrix output variable
  cout << Time() << ": Normalizing..." << endl;
  Normalize( norms, normed_matrix );
  
  
  // Find isolated contigs.  This is a way of handling sparse input data.
  // When the LinkMatrix is interpreted as an adjacency graph, not all nodes/contigs will be part of large components.  The extreme form of this is singleton
  // nodes: contigs with no adjacencies whatsoever, and therefore no links in the graph.  More generally, very small components (say, doubletons) almost
  // certainly represent noise at worst, and un-interpretable signal at best, so we filter them out.
  cout << Time() << ": FindIsolatedNodes" << endl;
  normed_matrix.FindIsolatedNodes( max_isolated_component_size );
  
  // Report on isolated contigs.
  int N_isolated = 0, N_non_isolated = 0, isolated_norm = 0, non_isolated_norm = 0;
  for ( int i = 0; i < size(); i++ )
    if ( normed_matrix._isolated_contigs[i] ) {
      N_isolated++;
      isolated_norm += norms[i];
    }
    else {
      N_non_isolated++;
      non_isolated_norm += norms[i];
    }
  cout << Time() << ": 'Isolated contigs' are in components of size <= " << max_isolated_component_size << " contigs.  N isolated = " << N_isolated << "; # RE sites = " << isolated_norm << ".  N non-isolated = " << N_non_isolated << "; # RE sites = " << non_isolated_norm << "." << endl;
  
  
  //normed_matrix.Print();
  
  /*
  // Find the coverage of each contig in the draft genome by shotgun reads, and use these coverages to re-weight the link densities.
  // The principle here is that contigs with wildly different coverage rates are unlikely to be from the same species.  (This doesn't take into account the
  // possibility of multi-species contigs, but true multi-species contigs will tend to have more coverage than all single-species contigs and will thus have
  // all of their true Hi-C links down-weighted.)
  vector<double> shotgun_covs = FindShotgunCoverages( shotgun_BAM, N );
  
  // Find the minimum non-zero coverage.  Use this in place of all zero coverages, to avoid dividing by zero.
  double min_shotgun_cov = INT_MAX;
  for ( size_t i = 0; i < N; i++ )
    if ( shotgun_covs[i] > 0 && shotgun_covs[i] < min_shotgun_cov )
      min_shotgun_cov = shotgun_covs[i];
  int N_changed = 0;
  for ( size_t i = 0; i < N; i++ )
    if ( shotgun_covs[i] == 0 ) {
      shotgun_covs[i] = min_shotgun_cov;
      N_changed++;
    }
  
  const int RATIO_POWER = 0; // TEMP: 0, 1, 2, 4 (0 does nothing)
  PRINT( RATIO_POWER );
  
  // Divide each link's density by the ratio of the shotgun read coverages of the two contigs.
  for ( LinkMatrixDouble::const_iterator1 it1 = normed_matrix.begin1(); it1 != normed_matrix.end1(); ++it1 )
    for ( LinkMatrixDouble::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2 ) {
      if ( *it2 == 0 ) continue;
      int i = it2.index1(), j = it2.index2();
      double ratio = shotgun_covs[i] / shotgun_covs[j];
      if ( ratio < 1 ) ratio = 1 / ratio;
      //PRINT3( i, j, ratio );
      normed_matrix(i,j) = normed_matrix(i,j) / pow( ratio, RATIO_POWER );
    }
  */
  
  
  
  // Apply the LimitNeighbors function from Jarvis-Patrick clustering.  The idea is to break links in the graph between contigs if the two contigs are not
  // in each other's list of top K nearest neighbors.  This removes outlier contigs in a local-density-independent manner, and also simplifies the graph so
  // that the subsequent hierarchical clustering is far faster and more straightforward.
  cout << Time() << ": LimitNeighbors" << endl;
  assert( normed_matrix._isolated_contigs.size() > 0 );
  normed_matrix.LimitNeighbors( K, true );
  assert( normed_matrix._isolated_contigs.size() > 0 );
  
}


// IntraClusterEnrichment: Input a vector indicating which cluster each contig belongs in.  Find the enrichment of links within a cluster.
// Skip clusters with total norm < min_cluster_norm.
// Another version of this same calculation is in MetaAssembly::FindIntraClusterEnrichment(), which inputs the clustering as a processed ClusteringResult.
template<class T> double
LinkMatrix<T>::IntraClusterEnrichment( const vector<int> & contig_clusterID, const vector<int> & contig_norms, const vector<int> & cluster_norms, const int min_cluster_norm ) const
{
  int64_t N_intra_cluster_links = 0, N_inter_cluster_links = 0;
  int64_t intra_cluster_len2 = 0, inter_cluster_len2 = 0;
  
  // Loop over all pairs of contigs that are clustered (into sufficiently large clusters.)
  // This isn't computationally efficient.
  for ( int i = 0; i < size(); i++ ) {
    int C1 = contig_clusterID[i];
    if ( C1 == -1 ) continue;
    if ( cluster_norms[C1] < min_cluster_norm ) continue;
    
    for ( int j = i+1; j < size(); j++ ) {
      int C2 = contig_clusterID[j];
      if ( C2 == -1 ) continue;
      if ( cluster_norms[C2] < min_cluster_norm ) continue;
      
      int len2 = contig_norms[i] * contig_norms[j];
      //int len2 = ( contig_norms[i] - 1 ) * ( contig_norms[j] - 1 ); // here, remove the +1 that's normally added to contig RE sites
      
      ( C1 == C2 ? N_intra_cluster_links : N_inter_cluster_links ) += (*this)(i,j);
      ( C1 == C2 ?   intra_cluster_len2  :   inter_cluster_len2  ) += len2;
    }
  }
  
  double intra_density = double( N_intra_cluster_links ) / double( intra_cluster_len2 );
  double inter_density = double( N_inter_cluster_links ) / double( inter_cluster_len2 );
  
  double enrichment = intra_density / inter_density;
  return enrichment;
}


// PrintMatrix: Print a visual representation of this matrix to stdout.  If it's a large matrix ( size >= 1000 ) then don't do anything.
template<class T> void
LinkMatrix<T>::Print() const
{
  if ( size() >= 1000 ) return; // shouldn't call this function on large matrices
  
  cout << "MATRIX:" << endl;
  
  for ( size_t i = 0; i < size(); i++ ) {
    for ( size_t j = 0; j < size(); j++ ) {
      cout << '\t';
      if ( (*this)(i,j) != 0 ) cout << (*this)(i,j);
      else if ( i <= j ) cout << '.';
    }
    cout << endl;
  }
}



// PrintSparse: Just print the locations where the sparse matrix is nonzero.
template<class T> void
LinkMatrix<T>::PrintSparse() const
{
  cout << "SPARSE MATRIX:" << endl;
  int N_values = 0;
  
  for ( typename LinkMatrix<T>::const_iterator1 it1 = this->begin1(); it1 != this->end1(); ++it1 )
    for ( typename LinkMatrix<T>::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2 )
      if ( *it2 != 0 ) {
	cout << "("<< it2.index1() << "," << it2.index2() << ") = " << *it2 << endl;
	N_values++;
      }
  
  cout << "SPARSE MATRIX TOTAL: " << size() << "x" << size() << ", N nonzero values = " << N_values << endl;
}





// WriteXGMML: Write this link matrix as an XGMML network file for visualization in Cytoscape.  Nodes are contigs, and edges are Hi-C links between contigs.
// If truth != NULL, give the nodes their true unique mappings.
// This is a lightweight format, with no additional information (node lengths, node names, cluster IDs, etc.)  Compare to MetaAssembly::WriteXGMML().
template<class T> void
LinkMatrix<T>::WriteXGMML( const string & XGMML_file, const T & MIN_LINK_WEIGHT, const TrueMapping * truth ) const
{
  cout << Time() << ": WriteXGMML!  Writing to file " << XGMML_file << " for Cytoscape to read..." << endl;
  
  int N_links_seen = 0, N_links_used = 0;
  
  // Open the XGMML file and write the XML header.
  ofstream XGMML( XGMML_file.c_str(), ios::out );
  
  XGMML << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<graph label=\"MetaAssembly of awesomeness\" directed=\"0\" xmlns=\"http://www.cs.rpi.edu/XGMML\">" << endl;
  
  vector<bool> seen( size(), false );
  vector< pair< pair<int,int>, T > > edges;
  
  // Loop over edges and determine which edges have sufficient link weight to be written.  Don't write them just yet.
  for ( int i = 0; i < size(); i++ ) {
    for ( int j = i+1; j < size(); j++ ) { // note that j > i
      
      T link_weight = (*this)(i,j);
      if ( link_weight == 0 ) continue;
      N_links_seen++;
      if ( link_weight < MIN_LINK_WEIGHT ) continue;
      N_links_used++;
      
      seen[i] = true;
      seen[j] = true;
      
      edges.push_back( make_pair( make_pair(i,j), link_weight ) );
    }
  }
  
  assert( N_links_seen == N_links );
  
  // For each contig, make a node.
  int N_nodes_used = 0;
  for ( int i = 0; i < size(); i++ )
    if ( seen[i] ) {
      XGMML << "  <node id=\"" << i << "\" label=\"" << i << "\">" << endl;
      N_nodes_used++;
      if ( truth != NULL )
	XGMML << "    <att name=\"ref_ID\" type=\"integer\" value=\"" << truth->QRefIDOnly(i) << "\"/>" << endl;
      XGMML << "  </node>" << endl;
    }
  
  
  // Write edges.
  for ( size_t e = 0; e < edges.size(); e++ ) {
    int i = edges[e].first.first;
    int j = edges[e].first.second;
    T link_weight = edges[e].second;
    XGMML << "  <edge source=\"" << i << "\" target=\"" << j << "\">" << endl;
    XGMML << "    <att name=\"link_weight\" type=\"real\" value=\"" << link_weight << "\"/>" << endl;
    XGMML << "  </edge>" << endl;
  }
  
  
  XGMML << "</graph>" << endl;
  XGMML.close();
  
  
  
  // Write a summary of what just happened.
  cout << Time() << ": WriteXGMML: Done!  Wrote " << N_nodes_used << " of " << size() << " nodes, and "
       << N_links_used << " of " << N_links << " edges (link weight >= " << MIN_LINK_WEIGHT << ")." << endl;
  
}
