#!/bin/bash

#SBATCH --partition rosa.p
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=0-10:00
#SBATCH --output=/user/rego3475/master_output/logs/6_make_trees.%j.out
#SBATCH --error=/user/rego3475/master_output/logs/6_make_trees.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ronja.roesner@uol.de


startdate=$(date '+%Y_%m_%d-%H_%M_%S')

mkdir $WORK/wd-tree_generation-$startdate
cd $WORK/wd-tree_generation-$startdate

mkdir treefiles-original treefiles-renamed treefiles-final logs trash
mkdir logs/automatic
mkdir logs/manual


# makes a phylotree for each alignment and renames it to be better readable
for file in master_input/aligned_hits;do
    locus_id=$(echo "${file%%-*}")
    echo === generating tree for $locus_id ===
    # creates a tree for this aligned hitfile
    iqtree2 -s $file -T 6 --tbe --alrt 10000    

    # renames tree branches to simplify analysis
    python3 ~/genotree/7_rename_trees.py -t $file.treefile -x ~/master_input/genome_master_library.xlsx

done

mv *aligned.fasta.treefile ./treefiles-original
mv *-renamed.treefile ./treefiles-renamed
mv *.fasta.* ./trash

# make a combined tree out of individual gene trees and run astral on it (program installed locally)
cat treefiles-renamed/*.treefile > treefiles-final/all-loci_combined.treefile
~/programs/ASTER-Linux/bin/astral-pro3 -t 8 -o treefiles-final/all-loci_astralpro3.treefile treefiles-final/all-loci_combined.treefile 2>astralpro3.log
~/programs/ASTER-Linux/bin/astral4 -t 8 -o treefiles-final/all-loci_astral4.treefile treefiles-final/all-loci_combined.treefile 2>astral4.log


mv astralpro3.log logs/automatic/astral3pro_$enddate.log
mv astral4.log logs/automatic/astral4_$enddate.log
mv genome_list.log logs/manual/genomelist_$enddate.log

enddate=$(date '+%Y_%m_%d-%H_%M_%S')
mkdir ~/master_output/phylo_trees/$enddate-trees
mv * ~/master_output/phylo_trees/$enddate-trees

cd ~/genotree
rm -r $WORK/wd-tree_generation-$startdate
