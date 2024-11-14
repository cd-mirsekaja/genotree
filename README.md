# genotree

Pipeline for making phylogenetic trees out of a dataset of genomes and exons. The scripts are written to run on a high performance cluster using SLURM.
The pipeline consists of multiple scripts working together.

0_make_alignments.sh - calls 1_process_genomes.sh, 3_align_loci.sh and 5_speciesinfo.py. Creates working directory and required subfolders while running.

1_process_genomes.sh - finds hits for each locus in specified genome file. Calls 2_find_hits.py

2_find_hits.py - turns nhmmer-tables into fasta files that contain hits for each locus

3_align_loci.sh - aligns hits for specified locus.

4_speciesinfo.py - gets general information for specified taxon

5_rate_alignments.sh - gives alignment ratings and removes bad ones

6_make_trees.sh - creates phylogenetic trees from the supplied alignments. Calls 7_rename_trees.py

7_rename_trees.py - renames tree tips to be more readable


