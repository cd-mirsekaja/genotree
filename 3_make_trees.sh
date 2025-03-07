#!/bin/bash

#SBATCH --partition rosa.p
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=0-24:00
#SBATCH --output=/user/rego3475/master_output/logs/3_make_trees.%j.out
#SBATCH --error=/user/rego3475/master_output/logs/3_make_trees.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ronja.roesner@uol.de

### Documentation ###
# This is a master script for generating phylogenetic supertrees
# from a set of multiple sequence alignments.
# Input:
# - a folder containing a number of FASTA files with MSAs
# 
# Output:
# - multiple folders containing individual gene trees
# - a folder containing consensus supertrees generated with ASTER


# load necessary modules
module load hpc-env/13.1
module load parallel/20230822-GCCcore-13.1.0
module load Python/3.11.3-GCCcore-13.1.0
module load IQ-TREE/2.2.2.7-gompi-2023a

# get the starting time
startdate=$(date '+%Y_%m_%d-%H_%M_%S')
echo === start date and time is $startdate ===

# make working directory and move into it
mkdir $WORK/wd-tree_generation-$startdate
cd $WORK/wd-tree_generation-$startdate

# make output directories
mkdir treefiles-original treefiles-renamed treefiles-final logs trash
mkdir logs/automatic

# set alignment file path and create genome list logfile
alignment_path=~/master_input/locus_CRY_full/aligned_filtered_realigned_files
touch genome_list.log

# set variable with alignment files
alignment_files=$(ls $alignment_path/*.fasta)
alignment_count=$(echo $filtered_alignments | wc -w)

# makes a phylotree for each alignment and renames it to be better readable
echo === generating individual gene trees at $(date '+%d.%m.%Y %H:%M:%S') ===
echo source for multiple sequence alignments $alignment_path

# create function for making the gene trees
generate_tree(){
    locus_id=$(echo "${1##*/}" | cut -d'-' -f1)
    echo === generating tree for $locus_id ===
    sbatch ~/genotree/3-1_make_gene_tree.sh $1
}

# run all tree generations in parallel
export -f generate_tree
echo "$alignment_files" | parallel --eta --jobs $alignment_count generate_tree

# waits for all trees to be generated by checking the squeue command for a certain term
while squeue -u $USER | grep -q "3-1_make"; do wait; done

# move iqtree output and renamed trees to their respective subfolders
mv $alignment_path/*.fasta.treefile ./treefiles-original
mv $alignment_path/*-renamed.treefile ./treefiles-renamed
mv $alignment_path/*.fasta.* ./trash

# make a combined tree out of individual gene trees and run astral on it (program installed locally)
echo === combining all gene trees and running ASTRAL at $(date '+%d.%m.%Y %H:%M:%S') ===
cat treefiles-renamed/*.treefile > treefiles-final/all-loci_combined.treefile
~/programs/ASTER-Linux/bin/astral-pro3 -t 8 -o treefiles-final/all-loci_astralpro3.treefile treefiles-final/all-loci_combined.treefile 2>astralpro3.log
~/programs/ASTER-Linux/bin/astral4 -t 8 -o treefiles-final/all-loci_astral4.treefile treefiles-final/all-loci_combined.treefile 2>astral4.log
~/programs/ASTER-Linux/bin/wastral -t 8 -o treefiles-final/all-loci_wastral.treefile treefiles-final/all-loci_combined.treefile 2>wastral.log
~/programs/ASTER-Linux/bin/caster-site -t 8 -o treefiles-final/all-loci_castersite.treefile treefiles-final/all-loci_combined.treefile 2>castersite.log
~/programs/ASTER-Linux/bin/caster-pair -t 8 -o treefiles-final/all-loci_casterpair.treefile treefiles-final/all-loci_combined.treefile 2>casterpair.log

# get the end time
enddate=$(date '+%Y_%m_%d-%H_%M_%S')

echo === finished ASTRAL. Moving files at $(date '+%d.%m.%Y %H:%M:%S') ===
# move log files to subfolder
mv astralpro3.log logs/automatic/astral3pro_$enddate.log
mv astral4.log logs/automatic/astral4_$enddate.log
mv wastral.log logs/automatic/wastral_$enddate.log
mv castersite.log logs/automatic/castersite_$enddate.log
mv casterpair.log logs/automatic/casterpair_$enddate.log
mv genome_list.log logs/genomelist_$enddate.log

# move all files to the output folder
mkdir ~/master_output/phylo_trees/$enddate-trees_CRY
mv * ~/master_output/phylo_trees/$enddate-trees_CRY

# remove the working directory
cd ~/genotree
#rm -r $WORK/wd-tree_generation-$startdate