// For documentation, see TrueMapping.h
#include "TrueMapping.h"
#include "TextFileParsers.h" // ParseTabDelimFile, TokenizeCSV, GetFastaNames


#include <assert.h>
#include <math.h>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm> // count

// Boost includes
#include <boost/dynamic_bitset.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

// Modules in ~/include (must add -L~/include and -lJ<module> to link)
#include "TimeMem.h"
#include "gtools/SAMStepper.h"




// TrueMapping constructor.
// Input a references CSV file describing a set of reference species, as well as a cache filename.
// This constructor uses blastn to align the assembly contigs onto the reference fastas.  It also caches the alignments to the cache file so that it doesn't
// have to do this work again in the future.  For species whose alignments are already stored in the cache file, just load from cache to save alignment time.
TrueMapping::TrueMapping( const string & assembly_fasta, const string & refs_csv_file, const string & refs_dir, const string & cache_file )
  : _query_names( GetFastaNames( assembly_fasta ) )
{
  _N_query_seqs = _query_names.size();
  
  _ref_names.clear();
  _ref_full_names.clear();
  _ref_fastas.clear();
  
  if ( !boost::filesystem::is_regular_file( refs_csv_file ) ) {
    cout << "ERROR: Can't find input CSV file `" << refs_csv_file << "'" << endl;
    exit(1);
  }
  
  // Parse the CSV file.  This gives us the number of reference species, the short and full names of each reference, and the location of the reference fasta.
  vector< vector<string> > csv_tokens;
  TokenizeCSV( refs_csv_file, csv_tokens );
  for ( size_t i = 0; i < csv_tokens.size(); i++ ) {
    vector<string> & tokens = csv_tokens[i];
    
    // Skip commented lines.
    assert( tokens.size() > 0 );
    assert( tokens[0].size() > 0 );
    if ( tokens[0][0] == '#' ) continue;
    
    
    // Do a bunch of sanity checks with verbose output, to make sure that any user error will immediately be caught.
    string error_str = "ERROR: Input CSV file `" + refs_csv_file + "', line " + boost::lexical_cast<string>(i);
    
    if ( tokens.size() != 3 ) {
      cout << error_str << ": Should see 3 tokens, instead saw " << tokens.size() << "." << endl;
      exit(1);
    }
    
    if ( find( _ref_names.begin(), _ref_names.end(), tokens[0] ) != _ref_names.end() ) {
      cout << error_str << ": Abbreviated name `" << tokens[0] << "' is already used by another line in this file.  These names should be unique!" << endl;
      exit(1);
    }
    
    if ( find( _ref_full_names.begin(), _ref_full_names.end(), tokens[1] ) != _ref_full_names.end() ) {
      cout << error_str << ": Species full name `" << tokens[1] << "' is already used by another line in this file.  These names should be unique!" << endl;
      exit(1);
    }
    
    string fasta = refs_dir + "/" + tokens[2];
    if ( !boost::filesystem::is_regular_file( fasta ) ) {
      cout << error_str << ": Can't find reference fasta file `" << fasta << "'." << endl;
      exit(1);
    }
    
    
    // Put the strings in the local structures.
    _ref_names     .push_back( tokens[0] );
    _ref_full_names.push_back( tokens[1] );
    _ref_fastas    .push_back( fasta );
    
  }
  
  _N_refs = _ref_names.size();
  
  cout << Time() << ": Creating a TrueMapping for an assembly on " << _N_refs << " reference ref_fastas" << endl;
  
  
  // Load the alignments!
  ReadAlignsAndWriteCache( assembly_fasta, cache_file );
}




void
TrueMapping::ReadAlignsAndWriteCache( const string & assembly_fasta, const string & cache_file )
{
  assert( _N_query_seqs > 0 );
  assert( _N_refs > 0 );
  
  // Initialize the alignment matrix, which indicates for each reference which of the query contigs align to it.  They are all initialized to 0 (false).
  _aligns.resize( _N_refs, boost::dynamic_bitset<>( _N_query_seqs, 0 ) );
  
  // Determine for each fasta whether or not it already has results in the cache file.
  vector<bool> has_cache( _N_refs, false );
  
  // If the cache file doesn't exist, this means there's no cache data, so we must align the assembly to all references.
  if ( !boost::filesystem::is_regular_file( cache_file ) ) {
    cout << Time() << ": Cache file " << cache_file << " does not exist.  Must align assembly to all references - this will take a while." << endl;
    for ( int i = 0; i < _N_refs; i++ )
      BlastAlign( assembly_fasta, i, cache_file );
    return;
  }
  
  
  // If the cache file does exist, read it and determine which fasta files have already had their alignments cached.
  
  // Make a lookup table of reference fasta file -> ref ID.
  map<string, int> ref_file_ID;
  for ( int i = 0; i < _N_refs; i++ )
    ref_file_ID[ _ref_fastas[i] ] = i;
  
  
  // List of which reference fastas have been seen in the cache file.
  vector<bool> ref_cached( _N_refs, 0 );
  
  // Read the cache file.
  vector< vector<string> > tokens;
  TokenizeFile( cache_file, tokens, true );
  
  // Each line in the file describes one reference genome, and should have exactly two tokens: [1] the name of the reference fasta file; and [2] a string of
  // _N_query_seqs integers (0/1) indicating whether or not each query contig aligned.
  for ( size_t i = 0; i < tokens.size(); i++ ) {
    assert( tokens[i].size() == 2 );
    map<string,int>::const_iterator it = ref_file_ID.find( tokens[i][0] );
    if ( it == ref_file_ID.end() ) continue; // the cache file may contain a reference fasta name not in _ref_fastas, which is fine but not useful
    int ref_ID = it->second;
    ref_cached[ref_ID] = true;
    
    // Read in the set of alignments from the string of 0's and 1's.  Note that inputting the dynamic_bitset directly from a string produces a reversed
    // result, with the lowest-indexed bits taken from the far right of the string.  This is fine, because that's how they were output in BlastAlign().
    _aligns[ref_ID] = boost::dynamic_bitset<>( tokens[i][1] );
    
    if ( _N_query_seqs != (int) _aligns[ref_ID].size() ) { // if this fails, the cache file describes multiple sets of queries (draft assemblies)
      cerr << "ERROR: The cache file at " << cache_file << " doesn't seem to have the right number of query sequences to match the assembly file at " << assembly_fasta << endl;
      exit(1);
    }
    assert( _N_query_seqs == (int) _aligns[ref_ID].size() );
  }
  
  // For all references that were not found in the cache, run BlastAlign.  This will find the alignments and also put them in the cache for the future.
  for ( int i = 0; i < _N_refs; i++ )
    if ( !ref_cached[i] )
      BlastAlign( assembly_fasta, i, cache_file );
  
}



// QOnRef: The most basic data function.  Given a query ID and a reference genome ID, return whether the query aligns.
bool
TrueMapping::QOnRef( const int qID, const int rID ) const
{
  assert( qID >= 0 && qID < _N_query_seqs );
  assert( rID >= 0 && rID < _N_refs );
  return _aligns[rID][qID];
}


// QOnRefOnly: Like QOnRef, but returns true only if this query aligns ONLY to this reference, and to no others.
bool
TrueMapping::QOnRefOnly( const int qID, const int rID ) const
{
  assert( qID >= 0 && qID < _N_query_seqs );
  assert( rID >= 0 && rID < _N_refs );
  for ( int i = 0; i < _N_refs; i++ )
    if ( ( i == rID ) != _aligns[i][qID] ) return false;
  return true;
}
  

// QRefIDOnly: If the query aligns to exactly one reference, returns that ref ID.  If it aligns to no refs, return -1; if multiple refs, return -2.
int
TrueMapping::QRefIDOnly( const int qID ) const
{
  assert( qID >= 0 && qID < _N_query_seqs );
  
  int rID = -1;
  for ( int i = 0; i < _N_refs; i++ )
    if ( _aligns[i][qID] ) {
      if ( rID == -1 ) rID = i;
      else return -2; // query aligns to multiple refs
    }
  
  return rID; // will either still be -1 (if qID doesn't align to any refs) or will be a ref ID >= 0
}


// QRefIDOnlyMY: Like QRefIDOnly, but for the MetaYeast scenario: considers refs #0,1,2,3 to be the same, returns 0 as a representative.
int
TrueMapping::QRefIDOnlyMY( const int qID ) const
{
  assert( qID >= 0 && qID < _N_query_seqs );
  
  int rID = -1;
  for ( int i = 0; i < _N_refs; i++ )
    if ( _aligns[i][qID] ) {
      if ( i < 4 ) rID = 0; // first four iterations are all S. cerevisiae species, so just mark this as aligning to SC-FY
      else if ( rID == -1 ) rID = i;
      else return -2; // query aligns to multiple refs, and not just multiple S. cerevisiae refs
    }
  
  return rID; // will either still be -1 (if qID doesn't align to any refs) or will be a ref ID >= 0
}



// QOnAnyRef: Return true iff the query aligns to at least one reference.
bool
TrueMapping::QOnAnyRef( const int qID ) const
{
  for ( int i = 0; i < _N_refs; i++ )
    if ( _aligns[i][qID] ) 
      return true;
  
  return false;
}


// RefsBitcode: Given a query ID, return an integer indicating which reference genome(s) it aligns to.  The integer is a set of consecutive bits: the i'th
// bit indicates whether the query aligns to reference #i (both bit and ref are 0-indexed.)  Hence the return value is in the range [0, 2^N_refs ).
int64_t
TrueMapping::RefsBitcode( const int qID ) const
{
  assert( qID >= 0 && qID < _N_query_seqs );
  assert( _N_refs < 8 * (int) sizeof( int64_t ) ); // make sure there aren't too many refs to prevent a return value (if this ever causes trouble, this function will need to be reimagined)
  
  int64_t code = 0;
  for ( int i = _N_refs - 1; i >= 0; i-- ) {
    code <<= 1;
    //if ( _aligns[i][qID] ) cout << i << '\t';
    if ( _aligns[i][qID] ) code++;
  }
  
  //cout << endl; PRINT2( qID, code ); cout << endl;
  return code;
} 



// ReportTargetCoverages: Input a SAM/BAM file with (shotgun) read alignments to the draft genome (the QUERY assembly).
// Report the average coverage of each query contig, and also report the average coverage of each reference genome.
void
TrueMapping::ReportTargetCoverages( const string & SAM_file ) const
{
  
  cout << Time() << ": Getting target lengths, coverages from file " << SAM_file << endl;
  vector<int>    contig_lens = TargetLengths( SAM_file );
  vector<double> contig_covs = TargetCoverages( SAM_file );
  
  assert( _N_query_seqs == (int) contig_lens.size() ); // if this fails, the SAM file's target sequences are not the same as this TrueMapping's query sequences
  assert( _N_query_seqs == (int) contig_covs.size() );
  
  cout << Time() << ": Reporting contig coverages" << endl;
  
  vector<int> ref_length( _N_refs, 0 ); // length of all contigs aligning to this ref
  vector<double> ref_cov_x_length( _N_refs, 0 ); // coverage * length of all contigs aligning to this ref
  
  // Report each contig's coverage, and which reference(s) it aligns to.
  // Also calculate, for each reference, the length-weighted average read coverage on all contigs that align to it.
  
  for ( int i = 0; i < _N_query_seqs; i++ ) {
    cout << "Contig #" << i << ":\t" << contig_covs[i];
    for ( int j = 0; j < _N_refs; j++ ) {
      bool on = _aligns[j][i];
      cout << '\t' << on;
      if ( on ) {
	ref_length[j] += contig_lens[i];
	ref_cov_x_length[j] += contig_covs[i] * contig_lens[i];
      }
    }
    cout << endl;
  }
  
  
  cout << Time() << ": Report" << endl;
  for ( int j = 0; j < _N_refs; j++ ) {
    double avg_cov = ref_cov_x_length[j] / ref_length[j];
    cout << "Species:\t" << _ref_names[j] << "\tShotgun coverage:\t" << avg_cov << endl;
  }
  
  cout << endl << endl;
}





// BlastAlign: Align the draft assembly to a reference fasta file.  Determine which contigs in the draft assembly align to the reference.
// This is time-consuming, especially for large reference genomes.  So once we have the results, we write them out to the cache file.
// This function requires the command-line programs `blastn` and `makeblastdb` (if necessary) to be available via the system() function.
// It also contains hard-wired parameters for blastn that affect alignment stringency and thus the interpretation of results.
void
TrueMapping::BlastAlign( const string & assembly_fasta, const int ref_ID, const string & cache_file )
{
  assert( ref_ID < _N_refs );
  const string ref_fasta = _ref_fastas[ref_ID];
  
  cout << Time() << ": BlastAlign: " << assembly_fasta << " onto " << ref_fasta << endl;
  
  
  // If necessary, use `makeblastdb` to create the BLAST database for the reference.
  string ref_blastdb = ref_fasta + ".blastdb";
  if ( !boost::filesystem::is_regular_file( ref_blastdb + ".nin" ) ||
       !boost::filesystem::is_regular_file( ref_blastdb + ".nhr" ) ||
       !boost::filesystem::is_regular_file( ref_blastdb + ".nsq" ) ) {
    string cmd = "makeblastdb -in " + ref_fasta + " -dbtype=nucl -out " + ref_blastdb;
    if ( system( cmd.c_str() ) != 0 ) {
      cout << "ERROR: makeblastdb did not complete successfully for some reason.  Can't make BLAST database for this reference, so can't align to it." << endl;
      exit(1);
    }
  }
  
  // Run blastn on this assembly.
  // TODO: make BLAST alignments more stringent so that the SC strains and SP will separate out
  // TODO: also, in the future, use the Blast API to do this natively in C++: http://www.ncbi.nlm.nih.gov/books/NBK7152/
  string blast_parameters = "-perc_identity 95 -evalue 1e-30 -word_size 50";
  // TEMP: alterations for the poplar metagenome assembly
  //blast_parameters = "-perc_identity 90 -evalue 1e-10 -word_size 20";
  
  string blast_file = cache_file + "." + boost::lexical_cast<string>(ref_ID) + ".blast";
  string cmd = "blastn -query " + assembly_fasta + " -db " + ref_blastdb + " " + blast_parameters + " -out " + blast_file + " -outfmt \"6 qseqid\" -num_threads 16";
  cout << "Command line:\t" << cmd << endl;
  if ( system( cmd.c_str() ) != 0 ) {
    cout << "ERROR: blastn did not complete successfully for some reason.  Hence, can't create BLAST alignments for reference " << ref_fasta << endl;
    exit(1);
  }
  
  // Make a list of flags (currently set to all false) for which query sequences have been found aligned to target.
  assert( ref_ID < (int) _aligns.size() );
  _aligns[ref_ID].resize( _N_query_seqs, false );
  
  // Make a lookup table of query sequence name -> ID.
  map<string,int> query_name_ID;
  for ( int i = 0; i < _N_query_seqs; i++ )
    query_name_ID[ _query_names[i] ] = i;
  
  
  
  // Because of how we ran blastn (-outfmt "6 qseqid") the output file is merely a single-column file containing a list of the query IDs for each alignment.
  // This is perfect, because all we want to know is which query sequences aligned to the reference genome - nothing about where on the reference they aligned.
  vector<string> aligned_names = ParseTabDelimFile<string>( blast_file, 0 );
  
  for ( size_t i = 0; i < aligned_names.size(); i++ ) {
    map<string,int>::const_iterator it = query_name_ID.find( aligned_names[i] );
    assert( it != query_name_ID.end() ); // if this fails, the contig names in the files <assembly_fasta> and <asssembly_fasta>.names are mismatched
    _aligns[ref_ID].set( it->second );
  }
  
  
  // Count and report the number of aligned sequences.
  int N_aligned = _aligns[ref_ID].count();
  cout << Time() << ": " << N_aligned << " of " << _N_query_seqs << " query sequences align to reference " << _ref_names[ref_ID] << " (using BLAST parameters: " << blast_parameters << ")" << endl;
  
  
  // Convert the result to a string of ones and zeros.  Output this string to the cache file (using 'append' rather than 'out' so as to not overwrite the file.)
  // Note that outputting the dynamic_bitset directly via << produces a reversed result, with the lowest-indexed bits on the far right.
  // This is fine, because that's how they'll be input later.
  ofstream out( cache_file.c_str(), ios::app );
  out << ref_fasta << '\t' << _aligns[ref_ID] << endl;
  out.close();
  
  
}
