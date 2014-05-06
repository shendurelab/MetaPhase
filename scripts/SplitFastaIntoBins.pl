#!/usr/bin/perl -w
use strict;

# Input a fasta file.  Split each contig into bins of length $N_BP_PER_BIN.  Useful for creating simulated reference genomes.
# This is NOT the same as the ChunkFasta.pl that splits a single assembly.fasta file into several files (for BLAST aligning) based on the contigs' lengths.


unless ( scalar @ARGV == 3 ) {
    print "\nSplitFastaIntoBins.pl: Split each contig in a fasta file into bins.  Useful for creating simulated reference genomes.\n\nSyntax: $0 <input-fasta> <output-fasta> <bin-size>\n\n";
    exit(1);
}
my ( $infile, $outfile, $BIN_SIZE ) = @ARGV;



# Constant parameters.
#my $outfile = "assembly.fasta";
my $N_BP_PER_LINE = 70; # number of bp in a line of a fasta file

# Human reference stuff.
my $use_canonical_human = 0; # set this to 1 if doing the human reference genome

#my $infile = "$ENV{'HOME'}/vol10/hg19/all/Homo_sapiens_assembly19.fasta";
my @canonical = qw / 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y/; # canonical chromosome names in human reference
my %canonical;
map {$canonical{$_}++} @canonical;





# SUBROUTINE
sub write_chrom_in_bins($$$)
{
    my ($chrom,$chrom_name,$fh) = @_;
    return unless $chrom_name; # skip the very first call, in which no chromosome has been read yet
    if ( $use_canonical_human ) { return unless $canonical{$chrom_name} }; # skip non-canonical chromosomes
    
    # Determine the number of bins in this chromosome.
    my $chrom_len = length $chrom;
    my $N_bins = int( $chrom_len / $BIN_SIZE ); # actually 1 less than the number of bins

    print localtime() . ": CHROMOSOME $chrom_name\tLength = $chrom_len\tN bins = $N_bins\n";
    
    # For each bin, print out a contig to the output file.
    foreach my $i (0..$N_bins) {
	print $fh ">${chrom_name}_bin$i\n"; # header line for this contig
	
	my $start = $i * $BIN_SIZE;
	my $contig = substr( $chrom, $start, $BIN_SIZE );
	my $contig_len = length $contig;
	
	# Split the contig into lines.
	my $N_lines = int( ($contig_len + $N_BP_PER_LINE - 1) / $N_BP_PER_LINE );
	foreach my $j (0..$N_lines-1) {
	    my $offset = $j * $N_BP_PER_LINE;
	    print $fh substr( $contig, $offset, $N_BP_PER_LINE ), "\n";
	}
    }
}




# MAIN

print "Bin size: $BIN_SIZE\n";
print "N bp per line in output file: $N_BP_PER_LINE\n";

# Open the fasta for reading.
my $bin_ID = 1;
die "$0: Can't find input file $infile: $!\n" unless -e $infile;
print localtime() . ": Opening input file $infile, writing to $outfile\n";
open IN, '<', $infile or die;
open OUT, '>', $outfile or die;


# Loop over all lines in the input fasta.
my $chrom_name = '';
my $chrom = '';

while (<IN>) {
    
    # If we're starting a new chromosome, then chop up the previous chrom into bins and write it to output.
    if ( /^>([\w\|]+)/ ) {
	&write_chrom_in_bins( $chrom, $chrom_name, \*OUT );
	$chrom = '';
	$chrom_name = $1;
    }
    
    # Otherwise, just append to this chromosome.
    else {
	chomp;
	$chrom .= $_;
    }
}


# Write the final chromosome.
&write_chrom_in_bins( $chrom, $chrom_name, \*OUT );

close IN or die;
close OUT or die;




