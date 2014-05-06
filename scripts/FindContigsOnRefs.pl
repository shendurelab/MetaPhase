#!/usr/bin/perl -w
use strict;

# FindContigsOnRefs.pl
#
# Determine which contigs in the MetaYeast draft assembly belong on which chromosomes, and tabulate and determine topline numbers for reporting.
#
# This script uses:
# -- Value of K used in making assembly.
# -- Directory containing BLAST output files at assembly.*.blast, which are produced by blast.sh.
# -- Reference strains' FastaSize files (hard-wired)
# -- Contig lengths.  These are taken from the FastaSize file.
#    NOTE: If the assembly was made with MetaVelvet, the lengths can be parsed directly out of the contigs' names so no additional files are needed.
#
# This script gets called by blast.sh for easy access.
#
# Josh Burton
# November 2013






# Command-line arguments.
unless ( @ARGV == 3 ) {
    print "Syntax: $0 <K> <out-dir-with-blast-files> <FastaSize-file>\n";
    exit;
}

my ( $K, $BLAST_dir, $FastaSize_file ) = @ARGV;

my $SPEEDUP = 0; # TEMP: option to speed up and get to the contigs, by skipping the slow reference coverage calculation

# Reference strains.
my @refs = ( 'SC-FY', 'SC-CEN', 'SC-RM', 'SC-SK1', 'SP', 'SM', 'SK', 'SB', 'NC', 'LW', 'LK', 'KL', 'KW', 'EG', 'SS', 'KP' );
my %ref_full_names = ( 'SC-FY' => 'S.cerevisiae-FY/S288C', 'SC-CEN' => 'S.cerevisiae-CEN.PK', 'SC-RM' => 'S.cerevisiae-RM11-1A', 'SC-SK1' => 'S.cerevisiae-SK1', 'SP' => 'S.paradoxus', 'SM' => 'S.mikatae', 'SK' => 'S.kudriavzevii', 'SB' => 'S.bayanus', 'NC' => 'N.castellii', 'LW' => 'L.waltii', 'LK' => 'L.kluyveri', 'KL' => 'K.lactis', 'KW' => 'K.wickerhamii', 'EG' => 'E.gossypii', 'SS' => 'S.stipitis', 'KP' => 'K.pastoris' );





sub N50(@) {
    
    my $total = 0;
    map { $total += $_ } @_;
    $total /= 2; # sum now has to be greater than this
    
    my $sum = 0;
    foreach (sort {$a <=> $b} @_) {
	$sum += $_;
	return $_ if $sum > $total;
    }
}


# contig_len_from_Velvet_name: Parse the name of a contig in a Velvet or MetaVelvet assembly to determine the contig's length.
sub contig_len_from_Velvet_name($) {
    die; # this is deprecated. use %contig_len instead.
    my @contig_name = split /_/, $_[0];
    #die "Weird contig name: $_[0]" unless scalar @contig_name == 6;
    my $contig_len = $contig_name[3] + ($K-1); # adjust the contig length by (K-1) - this is the wrong adjustment!
}



# Create a data structure that allows us to record the amount of coverage on each base in each chromosome of each reference genome.
# This requires reading in the reference genomes' FastaSize files.
print localtime() . ": Reading in reference FastaSize files and creating coverage data structure... (~20s)\n";
my %cov_on_ref; # 2D hash: ref name -> chrom name -> an array indicating coverage of each base on the chromosome
my %total_len; # total length of each reference
my $meta_len; # metagenome's entire length

foreach my $ref (@refs) {
    #print localtime() . ": $ref\n";
    
    my $file = "../refs/$ref/$ref.fasta.FastaSize";
    open IN, '<', $file or die "Can't find file $file: $!";
    
    # Parse each line of the FastaSize file and get the chromosome name and length.
    while (<IN>) {
	
	next if /TOTAL/; # skip last line
	die unless /^\s+(\d+)\s+(\S+)\n$/;
	my ($len,$chrom) = ($1,$2);
	$total_len{$ref} += $len;
	$meta_len += $len;
	
	next if $SPEEDUP;
	
	# Make the empty array for the coverage hash.
	my @cov = (0) x $len;
	$cov_on_ref{$ref}{$chrom} = \@cov;
    }
    
    close IN;
}





# Read the assembly's FastaSize file and make a mapping of contig name -> length.
my %contig_lens;

print localtime(). ": Reading assembly FastaSize file at $FastaSize_file...\n";
open IN, '<', $FastaSize_file or die;

while (<IN>) {
    last if /TOTAL/; # skip last line
    die "FastaSize file $FastaSize_file seems corrupt?" unless /\b(\d+)\s+(\S+)\b/;
    $contig_lens{$2} = $1;
}

close IN;




my %contig_ref_N; # will be a 2-D hash indicating the number of times each contig aligns to each ref, according to BLAST.
my %contig_ref_align_pos; # same as %contig_ref_N, but with each alignment's position instead of just tallying alignments.


# Open each BLAST file in turn to parse it.
foreach my $ref (@refs) {
    my $file = "$BLAST_dir/$ref.blast";
    print localtime() . ": $file\n"; $!++;
    open IN, '<', $file or die "Can't find file $file: $!";
    
    
    my $contig = '';
    
    # Parse the BLAST file, looking for lines that describe contigs by name and give their number of hits.
    while (<IN>) {
	#print "LINE: $_";
	
	my @tokens = split;
	if ( $tokens[0] eq '#' && $tokens[1] eq 'Query:' ) { # commented lines describing contig name
	    $contig = $tokens[2];
	}
	elsif ( $tokens[0] eq $contig ) { # line describing an alignment
	    #print "@tokens\n";
	    
	    # Find the region on the assembly contig and on the reference genome covered by this alignment.
	    my ($chrom,$c_start,$c_stop,$start,$stop) = @tokens[1,6,7,8,9];
	    ($start,$stop) = ($stop,$start) if $start > $stop;
	    
	    # Adjust the start/end regions to be zero-indexed, not one-indexed like the BLAST output.  Reduce the alignment by length (K-1) because that's what
	    # is truly represented by a de Bruijn graph contig.
	    $start--;
	    $c_start--;
	    $stop -= $K; # -1 b/c zero-indexed; -(K-1) because of de Bruijn graph
	    $c_stop -= $K;
	    
	    next if $start >= $stop; # this happens if the original length is not more than K
	    $contig_ref_N{$contig}{$ref}++; # mark this contig as covering this reference
	    my @c_range = ($c_start,$c_stop);
	    push @{$contig_ref_align_pos{$contig}{$ref}}, \@c_range;
	    next if $SPEEDUP;
	    
	    # Mark each base in the reference that is covered by this alignment.
	    #my @cov_array = exists $cov_on_ref{$ref}{$chrom} ? @{$cov_on_ref{$ref}{$chrom}} : ();
	    
	    map { $cov_on_ref{$ref}{$chrom}[$_]++ } ( $start..$stop ); # TODO: this is slooooooooow: increases runtime from 1s to ~90s
	    
	    #$cov_on_ref{$ref}{$chrom} = \@cov_array;
	    #print "$chrom\t$start\t$stop\t@cov_array\n";
	}
    }
    
    close IN;
}






# Now take each contig and look at which reference(s) it aligns to.  We now convert alignments to binary: either a contig aligns, or it doesn't.
# Tally up the total number and length of aligning contigs for each reference.
print localtime() . ": Tabulating results for each contig\n";

my %N_contigs_on_ref;
my %contig_lens_on_ref;
my @contig_lens = ();

my   %N_unique_on_ref; # for each reference, the number and length of contigs that align *only* to that reference
my %len_unique_on_ref;

my @contigs = sort keys %contig_lens;

my $meta_ref_coverage = 0;
my $meta_unique_ref_coverage = 0;


foreach my $contig (@contigs) {
    my $contig_len = $contig_lens{$contig};
    push @contig_lens, $contig_len;
    #print "$contig\t->\t$contig_len\n";
    
    my @refs_w_contig = ();
    
    foreach my $ref (@refs) {
	
	# Check whether this contig aligns to this reference.  If not, skip it.
	next unless $contig_ref_N{$contig}{$ref};
	
	# Add to tallies!
	$N_contigs_on_ref {$ref}++;
	push @{$contig_lens_on_ref{$ref}}, $contig_len;
	
	push @refs_w_contig, $ref;
    }
    
    # If this contig aligns to this reference and not to any others, note this fact.
    if ( scalar @refs_w_contig == 1 ) {
	$N_unique_on_ref  {$refs_w_contig[0]}++;
	$len_unique_on_ref{$refs_w_contig[0]} += $contig_len;
	$meta_unique_ref_coverage += $contig_len;
    }
    
}

my   $N50_contig_len = N50( @contig_lens );
my $total_contig_len = 0; map { $total_contig_len += $_ } @contig_lens;



# Report on how the reference genomes are covered by contigs.

print localtime() . ": Report on references!\n\n";
print "Total N contigs = ", scalar keys %contig_ref_N, "\n";
print "Total contig length = $total_contig_len\n";
print "N50 contig length = $N50_contig_len\n";

foreach my $ref (@refs) {
    
    # Tabulate the total length of all the covered stuff.
    my $total_cov = 0;
    foreach my $chr ( keys %{$cov_on_ref{$ref}} ) {
	map { $total_cov++ if $_ } ( @{$cov_on_ref{$ref}{$chr}} );
    }
    
    $meta_ref_coverage += $total_cov;
    
    # Find the N50 and total contig length on this reference.
    my $N50_len_on_ref = N50( @{$contig_lens_on_ref{$ref}} );
    my $total_len_on_ref = 0; map { $total_len_on_ref += $_ } @{$contig_lens_on_ref{$ref}};
    
    my $pct_cov = 100.0 * $total_cov / $total_len{$ref};
    printf "%-26sN contigs = $N_contigs_on_ref{$ref}\tN50 = $N50_len_on_ref\tTotal len = $total_len_on_ref\tTotal coverage = %-9d(%.2f%% of reference)\t  N contigs aligning only to this ref: $N_unique_on_ref{$ref}\t(len = $len_unique_on_ref{$ref})\n",
    $ref_full_names{$ref}, $total_cov, $pct_cov;
}

my $pct_meta_cov = 100.0 * $meta_ref_coverage / $meta_len;
printf "\nSummary:\nTotal length on all refs = $meta_ref_coverage (%.2f%%)\nTotal coverage by uniquely aligning contigs = $meta_unique_ref_coverage\n\n",
    $pct_meta_cov;



# Report on the contigs, sorting them by which reference genome(s) they align to.
print localtime(). ": Report on contigs!\n";

my   %N_contigs_on_refset; # map: set of reference genomes -> number and length of all contigs that align to exactly this set
my %len_contigs_on_refset;
my $chimeric_contigs_N = 0; # number and length of chimeric contigs
my $chimeric_contigs_len = 0;

foreach my $contig ( @contigs ) {
    
    # Find which reference(s) this contig aligns to.
    my @refs = sort grep { $contig_ref_N{$contig}{$_} != 0 } keys %{$contig_ref_N{$contig}};
    
    # Convert the set of references to a 'ref set' string.
    my $refset = join '_', scalar @refs, @refs;
    $N_contigs_on_refset{$refset}++;
    my $len = $contig_lens{$contig};
    $len_contigs_on_refset{$refset} += $len;
    
    
    # Determine whether this is a "chimeric" contig.  A chimeric contig is defined as one containing two different stretches, with at least 1kbp of sequence
    # each, for which the set of genomes they align to are different.  This requires parsing all of the contig's alignments to all references.
    my $CHIMERIC_BP = 1000;
    
    my $ref_align_pos = $contig_ref_align_pos{$contig};
    next unless $ref_align_pos;
    next unless scalar keys %$ref_align_pos > 1; # if this contig only ever aligns to one ref, it can't be chimeric
    
    # Tally up all alignments in the form of contig coverage.  Store the coverage as a set of "bits", indicating for each position on the contig, which
    # reference(s) it aligns to.  The first ref in $refs gets the 1 bit, the next gets the 2 bit, and so on.
    my @contig_bits = (0) x $len;
    #print "contig aligns from $contig:\n";
    my $ref_bit = 1;
    
    foreach my $ref ( sort keys %$ref_align_pos ) {
	foreach my $align ( @{$ref_align_pos->{$ref}} ) {
	    #print "ALIGN: $ref\t@$align\n";
	    map { $contig_bits[$_] += $ref_bit unless $contig_bits[$_] & $ref_bit } ( $align->[0]..($align->[1]-1) );
	}
	$ref_bit *= 2;
    }
    
    # Find the amount of sequence that falls into each type of coverage category, determined by which ref(s) it aligns to.
    my %bit_to_cov;
    foreach my $pos (0..$len-1) {
	next unless $contig_bits[$pos];
	#print "POS = $pos\tBIT = $contig_bits[$pos]\n";
	$bit_to_cov{ $contig_bits[$pos] }++;
    }
    
    
    # Examine the unique coverages and determine whether this contig is chimeric.
    #map { print "BIT: $_\tCOV: $bit_to_cov{$_}\n" } sort keys %bit_to_cov;
    my $N_unique_refs = scalar grep { $bit_to_cov{$_} >= $CHIMERIC_BP } keys %bit_to_cov;
    if ( $N_unique_refs > 1 ) {
	#print "N_UNIQUE_REFS = $N_unique_refs\n" ; $|++;
	$chimeric_contigs_N++;
	$chimeric_contigs_len += $len;
    }
    #print "\n";
    
    
}

printf "N chimeric contigs = $chimeric_contigs_N\nTotal chimeric contig length = $chimeric_contigs_len\tavg length = %.1f\n\n", $chimeric_contigs_len / $chimeric_contigs_N;

foreach my $refset ( sort {$len_contigs_on_refset{$b} <=> $len_contigs_on_refset{$a} } keys %len_contigs_on_refset ) {
    printf "N contigs = $N_contigs_on_refset{$refset}   \tlen = %-11s\t%-40s\n", $len_contigs_on_refset{$refset}, $refset;
}



print localtime() . ": Done!\n";
