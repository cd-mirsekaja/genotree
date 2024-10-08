# genotree

Pipeline for making phylogenetic trees out of a dataset of genomes and exons. The scripts are written to run on a high performance cluster using SLURM.

0_job_main_script.sh - calls 1_genomic_script.sh, 3_align_tree_script.sh and 5_speciesinfo.py. Creates working directory and required subfolders while running.

1_genomic_script.sh - finds hits for each locus in specified genome file. Calls 2_hmmer2fasta_script.py

2_hmmer2fasta_script.py - turns nhmmer-tables into fasta files that contain hits for each locus

3_align_tree_script.sh - aligns hits and creates a newick string for specified locus. Calls 4_treenaming_script.py

4_treenaming_script.py - renames tree tips to be more readable

5_speciesinfo.py - gets general information for specified taxon
