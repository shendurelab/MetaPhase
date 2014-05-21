## Scripts

The `scripts/` directory contains several scripts that I used to prepare the datasets for the paper.  I've included them here as a way of documenting what I did for.  I haven't modified them at all, so they won't run as-is (because, for example, they've got directory names from my filesystem hard-coded into them.)

- `assemble.sh`: A wrapper script for [IDBA-UD](i.cs.hku.hk/~alse/hkubrg/projects/idba_ud/index.html), the de novo metagenome assembler.  Used to create the M-Y (MetaYeast) draft assembly.
- `SplitFastaIntoBins.pl`: A script that splits all the contigs in a fasta file into bins of a specified size.  Used to create the M-3D simulated assembly.
- `FindContigsOnRefs.pl`: A script to determine which contigs in a draft assembly map to which of a set of reference genomes, and report summary stats.  Used to create the supplementary figure indicating draft contig coverage of the reference genomes.
- `blast.sh`: A wrapper script to call blastn.  Used in FindContigsOnRefs.pl.
- `align.iter.sh`: An iterative aligner.  Trims reads, tries to align them with [BWA](http://bio-bwa.sourceforge.net/), then trims them again and repeats.
