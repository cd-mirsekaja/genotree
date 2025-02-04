# genotree

Pipeline for making phylogenetic trees out of a dataset of genomes and exons. The scripts are written to run on a high performance cluster using SLURM.
The pipeline consists of multiple scripts working together. Input files and folders are set inside the scripts.


**Main Pipeline:**

1_make_alignments.sh

    - calls 1-1_process_genomes.sh, 1-3_align_loci.sh and utilities/speciesinfo.py
    
    - Creates working directory and required subfolders while running.


1-1_process_genomes.sh

    - calls 1-2_find_hits.py
    
    - finds hits for each locus in specified genome file.


1-2_find_hits.py

    - uses nhmmer-tables to create fasta files that contain hits for each locus


1-3_align_loci.sh

    - aligns hits for specified locus.


2_process_alignments.sh

    - calls 2-1_rate_alignments.sh, 2-2_filter_alignments.py and utilities/total_score.py
    
    - generates AliGROOVE-Scores for all alignments and returns .fasta-files that are filtered to not contain any sequences with a score below the supplied threshold (0.35)


2-1_rate_alignments.sh

    - generates AliGROOVE matrices with similarity scores for all sequence alignments.


2-2_filter_alignments.py

    - filters the supplied sequence alignments for the score threshold


3_make_trees.sh

    - calls 3-1_rename_trees.py
    
    - creates individual phylogenetic trees from the supplied alignments and combines them into a supertree


3-1_rename_trees.py

    - renames tree tips to make postprocessing easier



**Utilities:**

speciesinfo.py - gets general information for specified taxon

total_score.py - calculates total mean and median scores for supplied AliGROOVE matrices

rename_taxa.sh - calls rename_script.sh for all supplied .fasta-files

rename_script.py - renames genes in supplied .fasta file for processing with AliGROOVE
