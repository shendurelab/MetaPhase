/**************************************************************************************************************************************************************
 *
 * LinkMatrix.h
 *
 *
 * This module contains two struct types: LinkMatrixInt and LinkMatrixDouble.  These are convenience classes, consisting of template instantiations for sparse
 * Boost matrices, with the built-in assumption that the matrices are square (N1 = N2) and upper triangular ( matrix(i,j) = 0 unless i < j.)
 *
 * These classes also have built-in functions for link matrix analysis and manipulation, as well as member variables N_links and N_pairs (equal to the number
 * of nonzero elements and the total amount of data, respectively.)
 *
 * The LinkMatrix describes an alignment of Hi-C links to a draft metagenome assembly.  Each row/column in the matrix is a contig.  It does not keep track of
 * where the Hi-C reads are aligned to the contigs, so it is not good for representing alignments to a finished reference genome.
 *
 *
 *
 * Josh Burton
 * February 2014
 *
 *************************************************************************************************************************************************************/


#ifndef _LINK_MATRIX__H
#define _LINK_MATRIX__H

#include <boost/numeric/ublas/matrix_sparse.hpp> // mapped_matrix

//#include "TrueMapping.h"

#include <vector>
#include <string>
using namespace std;



// This is a Boost uBlas mapped_matrix, which is the type of matrix we use to represent Hi-C links between contigs.
//struct LinkMatrixInt;
//struct LinkMatrixDouble;

template<class T> struct LinkMatrix;
typedef LinkMatrix<int64_t> LinkMatrixInt;
typedef LinkMatrix<double> LinkMatrixDouble;




template<class T>
struct LinkMatrix : public boost::numeric::ublas::mapped_matrix<T>
{
  LinkMatrix<T>() { Resize(0); }
  
  int size() const { return this->size1(); }
  int has_data() const { return N_links != 0; } // returns true if LoadFromSAM() has been called
  
  // Resize: Resize (and initialize!) a matrix.  Forces the matrix to be square.  Use this, not the native mapped_matrix::resize() function.
  void Resize( size_t N );
  
  // LoadFromSAM: Load the data from a SAM/BAM file into this LinkMatrix.  This can be called repeatedly on the same LinkMatrix, but the SAM files must match.
  void LoadFromSAM( const string & SAM_file, const bool verbose = true );
  
  // Normalize: Normalize the link density in a matrix by contig norms (lengths).  Note that output is a LinkMatrix<double>.
  void Normalize( const vector<int> & norms, LinkMatrix<double> & normed_matrix ) const;
  
  // Bootstrap the Hi-C link data.  That is, re-create the dataset of links by using random sampling _with replacement_ from the current dataset.
  void Bootstrap( const bool set_random_seed = true ); 
  
  // LimitNeighbors: Remove links from the matrix so that every row has a maximum of K nonzero elements in it.  Specifically, remove all links L(a,b) from
  // the graph unless L is one of the top K highest numbers in both its row and its column: that is, it's one of the top K hits for both nodes a and b.
  // If reweight = true, also set the weights of all (non-removed) links to indicate the strength of overlapping nearest neighbors.
  // This step is also called "K-nearest neighbor sparsification" and is an important pre-processing step in the Jarvis-Patrick and SNN clustering algorithms.
  // See http://www-users.cs.umn.edu/~kumar/papers/kdd02_snn_28.pdf
  void LimitNeighbors( const int K, const bool reweight = true );
  
  // FindIsolatedNodes: Fill _isolated_contigs, which indicating for each contig whether it's part of a component of size <= max_isolated_component_size.
  void FindIsolatedNodes( const int max_isolated_component_size );
  
  // IsolateNormalizeLimit: Use the results of FindIsolatedNodes(), Normalize(), and LimitNeighbors() in that order.
  // Disconnect all nodes that are part of a component of size <= max_isolated_component_size; then normalize by contig norms; then run LimitNeighbors with K.
  // shotgun_BAM: A BAM file of shotgun reads aligned to the draft assembly.  Used to adjust link weights, on the assumption that contigs with very different
  //              coverages are less likely to belong to the same species cluster.  Can be an empty string, in which case no adjustment is done.
  void IsolateNormalizeLimit( const vector<int> & norms, const int max_isolated_component_size, const int K, const string & shotgun_BAM, LinkMatrix<double> & normed_matrix ) const;
  
  // IntraClusterEnrichment: Input a vector indicating which cluster each contig belongs in.  Find the enrichment of links within a cluster.
  // Another version of this same calculation is in MetaAssembly::FindIntraClusterEnrichment(), which inputs the clustering as a processed ClusteringResult.
  double IntraClusterEnrichment( const vector<int> & contig_clusterID, const vector<int> & contig_norms, const vector<int> & cluster_norms, const int min_cluster_norm = 0 ) const;
  
  // Print: Print a visual representation of this matrix to stdout.  If it's a large matrix ( size >= 1000 ) then don't do anything.
  // PrintSparse: Just print the locations where the sparse matrix is nonzero.
  void Print() const;
  void PrintSparse() const;
  
  // WriteXGMML: Write this link matrix as an XGMML network file for visualization in Cytoscape.  Nodes are contigs, and edges are Hi-C links between contigs.
  void WriteXGMML( const string & XGMML_file, const T & MIN_LINK_WEIGHT ) const;
  
  
  // Additional data structures.
  vector<bool> _isolated_contigs; // output of FindIsolatedNodes()
  int64_t N_links; // total number of nonzero elements in matrix
  T N_pairs; // sum of all elements in matrix
};





#endif
