#!/bin/bash

#SBATCH --partition rosa.p
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=0-24:00
#SBATCH --output=/user/rego3475/master_output/logs/3_make_trees.%j.out
#SBATCH --error=/user/rego3475/master_output/logs/3_make_trees.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ronja.roesner@uol.de

# load necessary modules
module load hpc-env/13.1
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

# makes a phylotree for each alignment and renames it to be better readable
echo === generating individual gene trees at $(date '+%d.%m.%Y %H:%M:%S') ===
for file in ~/master_input/all_hits_aligned_filtered/0_35;do
    locus_id=$(echo "${file##*/}" | cut -d'-' -f1)
    echo === generating tree for $locus_id ===
    # creates a tree for this aligned hitfile
    iqtree2 -s $file -T 6 --tbe --alrt 10000    

    # renames tree branches to simplify analysis
    python3 ~/genotree/3-1_rename_trees.py -t $file.treefile -d ~/master_input/genotree_master_library.db

done

# move iqtree output and renamed trees to their respective subfolders
mv *aligned.fasta.treefile ./treefiles-original
mv *-renamed.treefile ./treefiles-renamed
mv *.fasta.* ./trash

# make a combined tree out of individual gene trees and run astral on it (program installed locally)
echo === combining all gene trees and running ASTRAL at $(date '+%d.%m.%Y %H:%M:%S') ===
cat treefiles-renamed/*.treefile > treefiles-final/all-loci_combined.treefile
~/programs/ASTER-Linux/bin/astral-pro3 -t 8 -o treefiles-final/all-loci_astralpro3.treefile treefiles-final/all-loci_combined.treefile 2>astralpro3.log
~/programs/ASTER-Linux/bin/astral4 -t 8 -o treefiles-final/all-loci_astral4.treefile treefiles-final/all-loci_combined.treefile 2>astral4.log
~/programs/ASTER-Linux/bin/wastral -t 8 -o treefiles-final/all-loci_wastral.treefile treefiles-final/all-loci_combined.treefile 2>wastral.log
~/programs/ASTER-Linux/bin/caster-site -t 8 -o treefiles-final/all-loci_castersite.treefile treefiles-final/all-loci_combined.treefile 2>castersite.log
~/programs/ASTER-Linux/bin/caster-pair -t 8 -o treefiles-final/all-loci_casterpair.treefile treefiles-final/all-loci_combined.treefile 2>casterpair.log

echo === finished ASTRAL. Moving files at $(date '+%d.%m.%Y %H:%M:%S') ===
# move log files to subfolder
mv astralpro3.log logs/automatic/astral3pro_$enddate.log
mv astral4.log logs/automatic/astral4_$enddate.log
mv wastral.log logs/automatic/wastral_$enddate.log
mv castersite.log logs/automatic/castersite_$enddate.log
mv casterpair.log logs/automatic/casterpair_$enddate.log
mv genome_list.log logs/genomelist_$enddate.log

# get the end time and move all files to the output folder
enddate=$(date '+%Y_%m_%d-%H_%M_%S')
mkdir ~/master_output/phylo_trees/$enddate-trees
mv * ~/master_output/phylo_trees/$enddate-trees

# remove the working directory
cd ~/genotree
rm -r $WORK/wd-tree_generation-$startdate
