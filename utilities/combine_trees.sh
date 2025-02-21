#!/bin/bash

#SBATCH --partition rosa.p
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=0-4:00
#SBATCH --output=/user/rego3475/master_output/logs/utils-combine_trees.%j.out
#SBATCH --error=/user/rego3475/master_output/logs/utils-combine_trees.%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ronja.roesner@uol.de

# combines renamed gene trees into super trees using different verions of ASTER
tree_dir=${1}
savemod=${2}

cd /dss/work/rego3475/wd_test
mkdir treefiles-final/$savemod
out_dir=treefiles-final/$savemod


cat $tree_dir/*.treefile > $out_dir/all-loci_combined.treefile
~/programs/ASTER-Linux/bin/astral-pro3 -t 8 -o $out_dir/all-loci_astralpro3.treefile $out_dir/all-loci_combined.treefile 2>$out_dir/astralpro3.log
~/programs/ASTER-Linux/bin/astral4 -t 8 -o $out_dir/all-loci_astral4.treefile $out_dir/all-loci_combined.treefile 2>$out_dir/astral4.log
~/programs/ASTER-Linux/bin/wastral -t 8 -o $out_dir/all-loci_wastral.treefile $out_dir/all-loci_combined.treefile 2>$out_dir/wastral.log
~/programs/ASTER-Linux/bin/caster-site -t 8 -o $out_dir/all-loci_castersite.treefile $out_dir/all-loci_combined.treefile 2>$out_dir/castersite.log
~/programs/ASTER-Linux/bin/caster-pair -t 8 -o $out_dir/all-loci_casterpair.treefile $out_dir/all-loci_combined.treefile 2>$out_dir/casterpair.log