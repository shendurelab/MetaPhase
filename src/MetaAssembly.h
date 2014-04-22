/**************************************************************************************************************************************************************
 *
 * MetaAssembly.h
 *
 *
 * MetaAssembly is the central object of the MetaPhase program.  The MetaAssembly object describes a metagenome assembly.
 *
 * A MetaAssembly is initially loaded as a set of contigs, with no known relations between one another.  Then LoadFromSAM() is called, inputting one or more
 * SAM/BAM files containing Hi-C read pairs that are aligned to the contigs.  Each inter-contig read pair is interpreted as a link between the two contigs.
 * A graph is created (_link_matrix) in which the nodes are contigs and the edges are link counts between the contigs.  This graph can be used
 * as input for a clustering algorithm (Cluster) creating a ClusteringResult object.  It can be written in XGMML network format to a file for reading into
 * Cytoscape (WriteXGMML) or files for reading into Lachesis (WriteClusterFasta, WriteLachesisFiles).
 *
 *
 * Josh Burton
 * December 2013
 *
 *************************************************************************************************************************************************************/


#ifndef _META_ASSEMBLY__H
#define _META_ASSEMBLY__H

#include "LinkMatrix.h" // LinkMatrixInt
#include "TrueMapping.h"
#include "ClusteringResult.h"

#include <string>
#include <vector>
#include <iostream>

using namespace std;




class MetaAssembly
{
 public:
  
  /* CONSTRUCTORS */
  
  // Construct a MetaAssembly object from an assembly fasta file.
  MetaAssembly( const string & assembly_fasta, const string & RE_site_seq );
  
  // Load a TrueMapping object.
  void LoadTrueMapping( TrueMapping * truth ) { _truth = truth; }
  
  // Load a SAM/BAM file and put its links into this MetaAssembly.  (TODO: load multiple files at once; figure out what data structure to put this into)
  void LoadFromSAM( const string & SAM_file );
  
  // Find the average link density within clusters.
  void FindClusterDensities( const bool use_truth = true ) const;
  
  // Bootstrap the Hi-C link data.  That is, re-create the dataset of links by using random sampling _with replacement_ from the current dataset.
  void Bootstrap( const bool set_random_seed = true ) { _link_matrix.Bootstrap( set_random_seed ); }
  
  // Prepare the link matrix for clustering, by normalization etc.  This uses _link_matrix and fills _normed_link_matrix.
  void PrepLinkMatrix( const int K, const string & shotgun_BAM = "" );
  
  // Main clustering algorithm!
  void Cluster( const int N_clusters, const int MIN_CLUSTER_NORM );
  
  void SetClusters( const ClusteringResult & clusters ) { _clusters = clusters; }
  
  // FindMultiSpeciesContigs: Find contigs that could plausibly belong to more than one cluster, implying that their sequences appear in multiple species.
  void FindMultiSpeciesContigs();
  
  /* OUTPUT METHODS - these all require Cluster() to have been run */
  
  // Non-truth-based cluster analysis.
  // Another version of this same calculation is in LinkMatrix::IntraClusterEnrichment(), which inputs the cluster IDs of each contig as a vector.
  void FindIntraClusterEnrichment( ostream & out = cout ) const;
  
  // Analyze the contigs that aren't clustered.  Gives a more thorough report if truth data is available (_truth != NULL).
  void ReportUnclusteredContigs( ostream & out = cout ) const;
  
  // DrawFigure2a: Make an igraph output file that depicts this graph and clustering result as an awesome image.
  void DrawFigure2a( const string & OUT_PNG ) const;
  
  // Write this link matrix as an XGMML network file for visualization in Cytoscape.  Nodes are contigs, and edges are Hi-C links between contigs.
  // Clustering and truth information about each contig are included as attributes, and can be used within Cytoscape to change node size, color, etc.
  void WriteXGMML( const string & XGMML_file, const int MIN_LINK_WEIGHT = 1 ) const;
  
  // Write a fasta file representing all the contigs in one cluster.  The cluster chosen is cluster #ID - or, if is_ref_ID = true, the single largest cluster
  // corresponding to ref ID #ID.  If ID == -1 and !is_ref_ID, writes all of the contigs not placed in any cluster.
  void WriteClusterFasta( const string & output_fasta, const int ID, const bool is_ref_ID ) const;
  
  // Write a set of files describing one cluster, to be input to Lachesis.  The cluster chosen is the single largest cluster corresponding to ref #ref_ID.
  // The files include: all.GLM
  void WriteLachesisFiles( const string & Lachesis_dir, const int ref_ID ) const;
  
  
 private:
  
  // Helper function for DrawFigure2a().
  void DrawGraphvizEdges( ostream & dot, const vector< pair<int,int> > & edges, const string & color, const vector<bool> & contigs_to_use, const double weight ) const;
  
  
  
  
  // Fasta filename for the input draft assembly.  Used in the constructor, and also in WriteClusterFasta().
  string _assembly_fasta;
  
  // TrueMapping object, describing the truth placements of each contig on each reference genome.  Remains NULL if there is no TrueMapping available.
  TrueMapping * _truth;
  
  // Information on the contigs in this assembly.
  int _N_contigs;
  vector<string> _contig_names;
  vector<int> _contig_lengths, _contig_RE_sites; // RE_sites = restriction enzyme sites
  
  // The data structure of Hi-C links.  For two contigs i,j in [0,N_contigs), with i < j, _matrix(i,j) is the number of Hi-C links between i and j.
  // Note that _matrix(i,j) = 0 unless i < j.
  LinkMatrixInt _link_matrix;
  LinkMatrixDouble _normed_link_matrix; // matrix of link data, normalized and otherwise prepped for hierarchical clustering - see PrepLinkMatrix()
  
  // Data source file names.
  string _RE_sites_file;
  vector<string> _SAM_files;
  
  
  
  // Clustering results.  Determined in Cluster().  This object is publicly accessible.
 public:
  ClusteringResult _clusters;
};




#endif
