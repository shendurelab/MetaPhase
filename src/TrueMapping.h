/**************************************************************************************************************************************************************
 *
 * TrueMapping.h
 *
 * TODO: update doc for MetaPhase; this doc here is for Lachesis
 *
 * A TrueMapping is a list of all the contigs in a draft assembly, with an indication of where they map onto all known reference assemblies.  If you are
 * assembling a genome for which a reference assembly already exists, you'll want to set up a TrueMapping object so you can evaluate Lachesis' results.
 * On the other hand, if you're assembling a genome truly de novo, you'll never need to instantiate a TrueMapping.
 *
 * The motivation for a TrueMapping is as follows.  The Lachesis algorithm, as implemented in the GenomeLinkMatrix and ChromLinkMatrix classes, is designed
 * to take a set of contigs and scaffold them together into one or more scaffolds using Hi-C links that have been mapped onto them.  The best way to evaluate
 * the success of Lachesis is to scaffold together a draft assembly for an organism that also has a completed reference assembly, and then compare the
 * Lachesis-created scaffold with the "true" mapping of the draft assembly's contigs.  The true mapping is stored in the TrueMapping class, which is simply a
 * data structure designed for this purpose.
 *
 * The input to the TrueMapping class is a text file of BLAST output representing the alignment of the draft assembly to the true reference.  To create this
 * file, run blastn with -outfmt 7.  It may be a computational hassle to get blastn to run with a reasonable amount of time and memory usage.  To speed it up,
 * restrict the space of allowable alignments (here's a set of options that was used successfully on human: -perc_identity 99 -evalue 100 -word_size 50).
 *
 * Alternatively, a TrueMapping class can be created for a non-de novo assembly that consists of the reference genome split up into "bins" of equal length.
 *
 * Note that the contigs in the draft assembly may technically be scaffolds themselves, if they are given as contigs with long strings of interior N's to
 * represent gaps.  This is how ALLPATHS reports its draft scaffolds.  They are still treated as contigs by Lachesis.
 *
 *
 * Josh Burton
 * February 2013
 *
 *************************************************************************************************************************************************************/


#ifndef _META_PHASE_TRUE_MAPPING__H
#define _META_PHASE_TRUE_MAPPING__H

#include <assert.h>
#include <string>
#include <vector>
#include <iostream>

// Boost libraries
#include <boost/dynamic_bitset.hpp>

using namespace std;



// TODO: I don't really need this full-blown alignment object, all I'm going to record is whether or not a query aligned to a reference!
// Alignment: A helper struct for TrueMapping.
// Each Alignment object represents a single alignment from a query (a draft assembly contig) to a reference genome.  If a query aligns to multiple references,
// there will be one Alignment object for each such alignment.  However, there will only be one Alignment for each query-reference pairing, representing the
// single best alignment (highest E-value) BLAST found from that query to that reference.
struct Alignment
{
  int qID; // "query ID": the ID of the contig in the draft assembly that is aligned
  int rID; // "reference ID": the ID of the reference genome to which the query is aligned
  int tID; // "target ID": the ID of the sequence in reference genome #rID to which the query is aligned
  int start, stop; // the aligned region on the target sequence is from [start,stop).  If start > stop, then the alignment is RC and is from [stop,start).
};




class TrueMapping
{
 public:
  
  /* CONSTRUCTORS */
  
  // Default constructor.  This is only invoked when another object that contains a TrueMapping object is instantiated.
  TrueMapping() {}
  
  // TrueMapping constructor.
  // Input a references CSV file describing a set of reference species, as well as a cache filename.
  // This constructor uses blastn to align the assembly contigs onto the reference fastas.  It also caches the alignments to the cache file so that it doesn't
  // have to do this work again in the future.  For species whose alignments are already stored in the cache file, just load from cache to save alignment time.
  TrueMapping( const string & assembly_fasta, const string & refs_csv_file, const string & refs_dir, const string & cache_file );
  
  
  // ReadAlignsAndWriteCache: Fill the vector _aligns with bits indicating the ability of all contigs from the draft assembly to align onto all references.
  // Read as many alignments from the cache file as possible; for references not found in the cache file, run BlastAlign() (time-consuming) to get alignments.
  void ReadAlignsAndWriteCache( const string & assembly_fasta, const string & cache_file );
  
  /* 'GET' FUNCTIONS */
  
  int NQuerySeqs() const { return _N_query_seqs; }
  int NRefs() const { return _N_refs; }
  string RefName    ( const int rID ) const { assert( rID >= 0 && rID < _N_refs ); return _ref_names     [rID]; }
  string RefFullName( const int rID ) const { assert( rID >= 0 && rID < _N_refs ); return _ref_full_names[rID]; }
  
  // QOnRef: The most basic data function.  Given a query ID and a reference genome ID, return whether the query aligns.
  bool QOnRef( const int qID, const int rID ) const;
  // QOnRefOnly: Like QOnRef, but returns true only if this query aligns ONLY to this reference, and to no others.
  bool QOnRefOnly( const int qID, const int rID ) const;
  // QRefIDOnly: If the query aligns to exactly one reference, returns that ref ID.  If it aligns to no refs, return -1; if multiple refs, return -2.
  int QRefIDOnly( const int qID ) const;
  // QRefIDOnlyMY: Like QRefIDOnly, but for the MetaYeast scenario: considers refs #0,1,2,3 to be the same, returns 0 as a representative.
  int QRefIDOnlyMY( const int qID ) const;
  // QOnAnyRef: Return true iff the query aligns to at least one reference.
  bool QOnAnyRef( const int qID ) const;
  
  // RefsBitcode: Given a query ID, return an integer indicating which reference genome(s) it aligns to.  The integer is a set of consecutive bits: the i'th
  // bit indicates whether the query aligns to reference #i (both bit and ref are 0-indexed.)  Hence the return value is in the range [0, 2^N_refs ).
  int64_t RefsBitcode( const int qID ) const;
  
  // ReportTargetCoverages: Input a SAM/BAM file with (shotgun) read alignments to the draft genome (the QUERY assembly).
  // Report the average coverage of each query contig, and also report the average coverage of each reference genome.
  void ReportTargetCoverages( const string & SAM_file ) const;
  
  
 private:
  
  
  // BlastAlign: Align the draft assembly to a reference fasta file.  Determine which contigs in the draft assembly align to the reference.
  // This is time-consuming, especially for large reference genomes.  So once we have the results, we write them out to the cache file.
  // This function requires the command-line programs `blastn` and `makeblastdb` (if necessary) to be available via the system() function.
  // It also contains hard-wired parameters for blastn that affect alignment stringency and thus the interpretation of results.
  void BlastAlign( const string & assembly_fasta, const int ref_ID, const string & cache_file );
  
  
  /* DATA */
  
  int _N_refs; // number of species for which we have reference genomes; this is the size of many subsequent vectors including _species, _target_names

  // The following three vectors describe the reference species.  They must all have length _N_refs, and the corresponding entries should match.
  vector<string> _ref_names; // abbreviated tag names of reference species; must be unique but should also be short (e.g., "SB" for S. bayanus is fine)
  vector<string> _ref_full_names; // full binomial names of reference species (e.g., "Saccharomyces bayanus").  Used for some output figures.
  vector<string> _ref_fastas; // fasta files for each reference
  
  
  // Lists of contig names.
  int _N_query_seqs;
  vector<string> _query_names; // query names (draft assembly)
  
  // Alignments.  For each reference (of which there are _N_refs), there is a bitset of length _N_query_seqs, indicating whether or not each query contig
  // aligned to that target.
  vector< boost::dynamic_bitset<> > _aligns;

  // The alignments themselves.  The Alignment object is a local struct indicating query ID, reference ID, target ID on reference, and start/stop on target.
  //vector<Alignment> _aligns;
  
  // Indexes into the alignments.  These mappings are redundant with the _aligns vector, but allow us to quickly answer the following questions:
  // For a given query (i.e., a contig in the draft assembly), which reference(s) does it align to, and where?
  //vector< vector<int> > _query_align_IDs;
  // For a given reference assembly, which query sequences align to it?
  //vector< vector<int> > _ref_align_IDs;
  
  
};




#endif
