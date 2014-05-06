# Iterative alignment!  Align the Hi-C reads to a draft assembly.
# We do this iteratively: trim the reads to 125bp; try to find perfect 125-mer alignments; trim the unaligned reads to 10bp and try again; then 75bp and 50bp.
# The input fastq files are made by trim.sh.
source $HOME/.bashrc


BWA_VERSION=0.6.2
BWA=/net/gs/vol3/software/modules-sw/bwa/$BWA_VERSION/Linux/RHEL6/x86_64/bin/bwa # DO NOT use the latest version of bwa - it's unnecessarily buggy

BAM2FASTQ=/net/shendure/vol10/jnburton/extern/bam2fastq-1.1.0/bam2fastq



# run_cmd: Run the command in $1.  If a $2 is given, pipe all output of $1 to file $2.
# Make sure to input $1 and $2 as quoted strings!
function run_cmd {
    echo `date`: $1
    #return # de-comment if this is a dry run
    if [ $# == 1 ]; then
	$1
    else
	$1 >& $2
    fi
}




# Scenario switch!
SCENARIO=M2
ALN_OPTS="-t 8 -n 0" # require a perfect match! (MY, poplar)
ALN_OPTS="-t 24" # generous alignment (M2)

# Read library to align.
LIB=M2.H3_pelleted # TEMP: Hind3, Nco1 (for MY); M2.H3, M2.N1, M2.H3_pelleted, M2.N1_pelleted (for M2)
RE_SITE=AAGCTT # HindIII = AAGCTT; NcoI = CCATGG


# Align the reads to the draft assembly. - MetaYeast or M2
#REF=../assembly/ASM2.fasta
#DIR=MY/to_ASM2

ASSEMBLY=sim.10k
REF=../assembly/$SCENARIO/$ASSEMBLY.fasta
DIR=$SCENARIO/to_$ASSEMBLY

# Align to an individual species' reference, or the meta-reference
#SPECIES=SS # species abbreviation (e.g., SC-FY), or 'combined_ref'
#REF=../refs/$SPECIES/$SPECIES.fasta
#DIR=MY/to_$SPECIES








# Set up the directory.
mkdir -p $DIR
ln -sf  ../$LIB.trim100.1.fq $DIR
ln -sf  ../$LIB.trim100.2.fq $DIR


echo `date`: start


for LEN in 100 75 50 ; do
    
    HEAD=$DIR/$LIB.trim$LEN
    
    # Run bwa aln and bwa sampe to align this pair of fastq files and create a SAM file.
    # The alignment requires a perfect match, so a lot of reads won't make it.  That's ok for now, we'll trim them soon.
    
    run_cmd "$BWA aln $ALN_OPTS $REF $HEAD.1.fq -f $HEAD.1.sai"
    run_cmd "$BWA aln $ALN_OPTS $REF $HEAD.2.fq -f $HEAD.2.sai"
    run_cmd "$BWA sampe -s $REF $HEAD.1.sai $HEAD.2.sai $HEAD.1.fq $HEAD.2.fq -f $HEAD.sam"
    
    # SAM -> BAM, flagstat
    run_cmd "make_bam_and_flagstat $HEAD.sam"
    # Grab out the read pairs in which both reads mapped.  These are the 'good' read pairs that we can use.
    # (Reads mapping with MQ=0 is a problem, but not one that can be solved by trimming the reads, which is what we're about to do with the unmapped pairs.)
    run_cmd "samtools view -bh -F 12 $HEAD.bam -o $HEAD.good.bam"
    run_cmd "samtools flagstat $HEAD.good.bam" "$HEAD.good.flagstat"
    
    if [ $LEN == 50 ] ; then break ; fi # end the last iteration here
    
    # Take the 'unaligned' read pairs (i.e., read pairs in which the reads didn't both align).  Convert them to a pair of fastq files, then trim them by 25 bp.
    NEXT_LEN=`expr $LEN - 25`
    NEXT_HEAD=$DIR/$LIB.trim$NEXT_LEN
    run_cmd "SamToFastqs.pl -P -t $NEXT_LEN $NEXT_HEAD.1.fq $NEXT_HEAD.2.fq $HEAD.sam"
done


# Merge the resulting set of 'good' bam files.
run_cmd "samtools merge -fn $DIR/$LIB.all_lens.bam $DIR/$LIB.trim100.good.bam $DIR/$LIB.trim75.good.bam $DIR/$LIB.trim50.good.bam"
run_cmd "samtools flagstat $DIR/$LIB.all_lens.bam" "$DIR/$LIB.all_lens.flagstat"
run_cmd "samtools view -h $DIR/$LIB.all_lens.bam -o $DIR/$LIB.all_lens.sam"

# REduce the data - that is, filter the BAM files to include only reads within 500bp of a RE site.
# This is also done in REduce.sh.
run_cmd "$HOME/L/git/make_bed_around_RE_site.pl $REF $RE_SITE 500"
run_cmd "bedtools intersect -abam $DIR/$LIB.all_lens.bam -b $REF.near_$RE_SITE.500.bed" "$DIR/$LIB.all_lens.REduced.bam"
run_cmd "samtools flagstat $DIR/$LIB.all_lens.REduced.bam" "$DIR/$LIB.all_lens.REduced.flagstat"

echo `date`: Done!
