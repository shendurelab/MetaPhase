#!/bin/bash -f

# Assemble the MetaYeast sequence library using IDBA-UD!
# This follows the example given in the IBDA README file at <extern>/ibda-1.1.1/README.
#
# Josh Burton
# January 2014

source ~/.bashrc


K=55
OUT_DIR=out

# Fragment libraries.
FQ1=../fastq/MetaYeast_1.fastq
FQ2=../fastq/MetaYeast_2.fastq
# Jumping libraries.
#JUMP_S1_FQ1=../jump_libs/S1.1.fastq # old jumping fastq's - short reads
#JUMP_S1_FQ2=../jump_libs/S1.2.fastq
#JUMP_S2_FQ1=../jump_libs/S2.1.fastq
#JUMP_S2_FQ2=../jump_libs/S2.2.fastq
#JUMP_S3_FQ1=../jump_libs/S3.1.fastq
#JUMP_S3_FQ2=../jump_libs/S3.2.fastq
JUMP_FQ1=../jump_libs/3kb_jump_combined_r1.fq
JUMP_FQ2=../jump_libs/3kb_jump_combined_r3.fq

IDBA_DIR=/net/shendure/vol10/jnburton/extern/idba-1.1.1/bin

# Convert the paired fastq files into interlaced fasta files.  Runtime: 13min for the fragment libraries; 1min total for the short-read jumping libraries.
#echo `date`: fq2fa
#$IDBA_DIR/fq2fa --merge $FQ1 $FQ2 frag.fasta
#echo `date`: fq2fa jump1
#$IDBA_DIR/fq2fa --merge $JUMP_FQ1 $JUMP_FQ2 jump.fasta
#echo `date`: fq2fa jump1
#$IDBA_DIR/fq2fa --merge $JUMP_S1_FQ1 $JUMP_S1_FQ2 jump1.fasta
#echo `date`: fq2fa jump2
#$IDBA_DIR/fq2fa --merge $JUMP_S2_FQ1 $JUMP_S2_FQ2 jump2.fasta
#echo `date`: fq2fa jump3
#$IDBA_DIR/fq2fa --merge $JUMP_S3_FQ1 $JUMP_S3_FQ2 jump3.fasta

# What needs to happen: more stringent assembly! - reduce chimerism, increase % uniquely aligning contigs
# Options to consider tweaking (w/defaults): add --pre_correction, --min_count 2, --min_support 1
# others: --mink 20, --maxk 100, --step 20, --inner_mink 10, --inner_step 5
OPTS="--pre_correction --maxk 60 --step 10"
OUT_DIR=$OUT_DIR.pc.20-10-60
echo `date`: idba_ud
#tm $IDBA_DIR/idba_ud $OPTS --read frag.fasta --read_level_2 jump1.fasta --read_level_3 jump2.fasta --read_level_4 jump3.fasta -o $OUT_DIR
tm $IDBA_DIR/idba_ud $OPTS --read frag.fasta --read_level_2 jump.fasta -o $OUT_DIR
echo `date`: done!
