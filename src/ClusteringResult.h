/**************************************************************************************************************************************************************
 *
 * ClusteringResult.h
 *
 *
 * A ClusteringResult is a data structure that holds the result of a clustering algorithm, in which contigs are assigned to clusters.
 * It is initially produced by the MetaAssembly::Cluster() function, and can be modified by MetaAssembly::FindMultiSpeciesContigs().
 * Note that any one contig may not be assigned to any clusters, or it may be assigned to multiple clusters if it is deemed a multi-species contig.
 *
 *
 * Josh Burton
 * January 2014
 *
 *************************************************************************************************************************************************************/


#ifndef _CLUSTERING_RESULT__H
#define _CLUSTERING_RESULT__H

#include "TrueMapping.h"

#include <string>
#include <vector>
#include <iostream>


using namespace std;





class ClusteringResult
{
 public:
  
  /* CONSTRUCTORS */
  
  ClusteringResult();
  // Construct a ClusteringResult with a set of contigs but no cluster assignments.
  ClusteringResult( const vector<int> & contig_lengths, const int N_clusters = 0 );
  // Construct a ClusteringResult object from the immediate result of clustering.
  ClusteringResult( const vector<int> & contig_to_clusterID, const vector<int> & contig_lengths );
  // Construct a ClusteringResult object from a ClusteringResult previously saved to file.
  ClusteringResult( const string & filename ) { ReadFile( filename ); }
  
  // ReadFile, WriteFile: Read/write a clustering result from file.
  // This alows you to evade the runtime cost of running the MetaAssembly::Cluster() function, once you've already run it with the same dataset and parameters.
  void ReadFile ( const string & filename );
  void WriteFile( const string & filename ) const;
  
  
  // Add a contig to a cluster.  Remove it from whatever cluster it was in before.  This function can be used to create a new cluster at #_N_clusters.
  void SetCluster( const int contig_ID, const int cluster_ID );
  
  // ReorderClustersByRefs: Renumber the cluster IDs in increasing order of PluralityRefID().
  // This uses truth information, but the cluster order doesn't really matter and is just for convenience/visual appeal, so this isn't really cheating.
  void ReorderClustersByRefs( const TrueMapping & truth );
  
  
  /* QUERY FUNCTIONS */
  
  bool empty() const { return ( _N_contigs == 0 ); } // will return true iff this object was created with the null constructor and hasn't been re-created
  int NContigs()  const { return _N_contigs; }
  int NClusters() const { return _N_clusters; }
  int cluster_ID( const int contig_ID )     const { assert(  contig_ID >= 0 &&  contig_ID <  _N_contigs ); return _contig_clusters[contig_ID]; }
  int ClusterSize( const int cluster_ID )   const; // ClusterSize: return the number of contigs in this cluster (currently computationally inefficient)
  int ClusterLength( const int cluster_ID ) const { assert( cluster_ID >= 0 && cluster_ID < _N_clusters ); return _cluster_lengths[cluster_ID]; }
  
  // PluralityRefID: Determine the reference ID to which a plurality of the contigs in this cluster align.
  int    PluralityRefID  ( const TrueMapping & truth, const int cluster_ID ) const;
  string PluralityRefName( const TrueMapping & truth, const int cluster_ID ) const;
  int LargestClusterWithPluralityRefID( const TrueMapping & truth, const int ref_ID ) const; // return -1 if no cluster has this plurality ref ID
  
  /* OUTPUT FUNCTIONS */
  
  // Print: Output summary statistics.
  void Print( const TrueMapping & truth, ostream & out = cout ) const;
  
  // Truth-based evaluation of the clusters in this ClusteringResult.  Produce the truth heatmap and chart files.  If plot = true, plot the heatmap image too.
  // If MY = true, assume this is a MetaYeast run, and hide SC-non-FY and KP.
  void DrawChart( const TrueMapping & truth, const bool require_unique_alignment, const bool plot, const bool MY, const int MIN_CLUSTER_LEN = 0 ) const;
  
  // DrawFigure2bc: Call DrawChart to create the chart and heatmap files, then call a script in ../figs to generate Figures 2b and 2c of the paper.
  void DrawFigure2bc( const TrueMapping & truth, const bool MY ) const;
  
  
 private:
  
  // FindClusterLengths: Input _contig_lengths and _contig_clusters; fill _cluster_lengths.  This is purely a bookkeeping operation.
  void FindClusterLengths();
  
  
  // Information on the contigs in this assembly.
  int _N_contigs;
  vector<int> _contig_lengths;
  
  
  // Information on the clusters.
  int _N_clusters;
  vector<int> _contig_clusters; // contig ID -> derived cluster ID (-1 if unclustered)
  vector<int> _cluster_lengths; // cluster ID -> length of all contigs in cluster
  
  
  
  // Friend function declarations (see below for documentation on these functions.)
  friend ClusteringResult MergeClusteringResults1( const ClusteringResult & cr1, const ClusteringResult & cr2 );
  friend ClusteringResult MergeClusteringResults2( const vector<ClusteringResult> & crs, const int new_N_clusters );
  friend ClusteringResult MergeClusteringResults3( const vector<ClusteringResult> & crs, const int new_N_clusters );
};





// MergeClusteringResults: Combine two sets of ClusteringResults into a single ClusteringResult by matching up the clusters and finding intelligent ways to
// split/merge/hem/haw little bits of cluster.
// These functions are designed as an intelligent way to merge ClusteringResults produced by Hi-C datasets with differing restriction sites.
// Method (MergeClusteringResults1): straightforward "intersection" method.
// Method (MergeClusteringResults2): Monti et al, http://www.cs.utexas.edu/users/ml/biodm/papers/MLJ-biodm2.pdf
// Method (MergeClusteringResults2): Dimitriadou et al, http://epub.wu.ac.at/94/1/document.pdf
ClusteringResult MergeClusteringResults1( const ClusteringResult & cr1, const ClusteringResult & cr2 );
ClusteringResult MergeClusteringResults2( const vector<ClusteringResult> & crs, const int new_N_clusters );
ClusteringResult MergeClusteringResults3( const vector<ClusteringResult> & crs, const int new_N_clusters );


#endif
