#!/bin/bash

# Align a set of contigs (a draft assembly in progress) to all 16 yeast reference genomes.
# BLAST (not cross_match, not BWA) is the best aligner for this.
#
# After aligning, use FindContigsOnRefs.pl to create a report about these alignments at FindContigsOnRefs.out.  This helps evaluate the quality of the de novo
# assembly and compare it to other assemblies.
#
# Josh Burton
# November 2013

#DIR=/net/shendure/vol10/jnburton/src/MetaYeast/assembly
DIR=.



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






# User-specified parameters.
K=63
ASSEMBLY_FASTA=$DIR/velvet/$K/contigs.fa
BLAST_OUT_DIR=$DIR/blast/velvet.K$K

# Try the IDBA assembly.
K=60

# Try the SPAdes assembly.
#K=55
#ASSEMBLY_FASTA=$DIR/SPAdes/out/scaffolds.fasta
#BLAST_OUT_DIR=$DIR/blast/SPAdes.scaffolds


# Different version of the IDBA assembly.
for CASE in contig scaffold ; do
#for CASE in contig scaffold scaffold-level-2 scaffold-level-3 scaffold-level-4 ; do
    
    ASSEMBLY_FASTA=$DIR/IDBA-UD/out.pc.20-10-$K/$CASE.fa
    BLAST_OUT_DIR=$DIR/blast/IDBA.pc.20-10-$K.$CASE
    
    mkdir -p $BLAST_OUT_DIR
    
    
    YEASTS="SC-FY SC-RM SC-CEN SC-SK1 SP SM SK SB NC LW LK KL KW EG SS KP"
    
    for YEAST in $YEASTS ; do
	REFDB=$DIR/../refs/$YEAST/$YEAST.fasta.blastdb
        # If you change the e-value here, remember to also change it in ~/MY/src/TrueMapping.cc!
	CMD="blastn -query $ASSEMBLY_FASTA -db $REFDB -perc_identity 95 -evalue 1e-30 -word_size 50 -out $BLAST_OUT_DIR/$YEAST.blast -outfmt 7"
	CMD="$CMD -num_threads 24"
	run_cmd "$CMD"
    done
    
    
    # Run FastaSize if necessary.
    if [ ! -e $ASSEMBLY_FASTA.FastaSize ] ; then
	echo `date`: FastaSize
	run_cmd "FastaSize $ASSEMBLY_FASTA" "$ASSEMBLY_FASTA.FastaSize"
    fi
    
    echo `date`: "FindContigsOnRefs.pl > $BLAST_OUT_DIR/FindContigsOnRefs.out (runtime: ~2m)"
    run_cmd "FindContigsOnRefs.pl $K $BLAST_OUT_DIR $ASSEMBLY_FASTA.FastaSize" "$BLAST_OUT_DIR/FindContigsOnRefs.out"
done


echo `date`: Done!
