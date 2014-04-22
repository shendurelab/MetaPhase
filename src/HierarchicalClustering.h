/**************************************************************************************************************************************************************
 *
 * HierarchicalClustering
 *
 * This module implements the agglomerative hierarchical clustering algorithm.
 * The main function is AgglomerativeHierarchicalClustering().  It's modeled after the function GenomeLinkMatrix::AHClustering() in the Lachesis code base.
 *
 * Josh Burton
 * January 2014
 *
 *************************************************************************************************************************************************************/

#ifndef _HIERARCHICAL_CLUSTERING__H
#define _HIERARCHICAL_CLUSTERING__H


#include "LinkMatrix.h" // LinkMatrixInt, LinkMatrixDouble

#include <boost/numeric/ublas/matrix.hpp> // matrix

#include <vector>
using namespace std;




/* AgglomerativeHierarchicalClustering: The main clustering algorithm used by MetaPhase.
 *
 * Inputs:
 *
 * matrix_orig: An NxN matrix indicating the amount of connectivity (e.g., Hi-C links) between two nodes i,j.  Note that the matrix is expected to be upper
 *              triangular: that is matrix(i,j) = 0 for i >= j.  This matrix is used only for reporting intra-cluster enrichment.
 * matrix_in:   A normalized form of matrix_orig.  This is the version used for the actual algorithm.
 * norms:       A vector of length N, giving the lengths of the contigs, either in bp or in number of restriction sites.
 * N_clusters:  The number of clusters to create, i.e., the number of distinct species to call
 * MIN_CLUSTER_NORM: clusters with a total norm smaller than this are not counted toward N_clusters and are discarded at the end
 *
 */
vector<int>
AgglomerativeHierarchicalClustering( const LinkMatrixInt & matrix_orig, const LinkMatrixDouble & matrix_in, const vector<int> & norms, const int N_clusters, const int MIN_CLUSTER_NORM );


// A simpler version of AgglomerativeHierarchicalClustering that just performs the basic algorithm.  Used in the MergeClusteringResults2() function.
vector<int>
AgglomerativeHierarchicalClusteringSimple( const boost::numeric::ublas::matrix<double> & matrix_in, const int N_clusters );




#endif
