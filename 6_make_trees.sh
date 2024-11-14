

echo === generating tree for $locus_id ===
# creates a tree for this aligned hitfile
iqtree2 -s $locus_id-aligned.fasta -T 6 --tbe --alrt 10000

# renames tree branches to simplify analysis
python3 ~/genotree/7_rename_trees.py -t $locus_id-aligned.fasta.treefile -x ~/master_input/genome_master_library.xlsx

# make a combined tree out of individual gene trees and run astral on it (program installed locally)
cat treefiles-renamed/*.treefile > treefiles-final/all-loci_combined.treefile
~/programs/ASTER-Linux/bin/astral-pro3 -t 8 -o treefiles-final/all-loci_astralpro3.treefile treefiles-final/all-loci_combined.treefile 2>astralpro3.log
~/programs/ASTER-Linux/bin/astral4 -t 8 -o treefiles-final/all-loci_astral4.treefile treefiles-final/all-loci_combined.treefile 2>astral4.log


mv astralpro3.log logs/automatic/astral3pro_$enddate.log
mv astral4.log logs/automatic/astral4_$enddate.log
mv genome_list.log logs/manual/genomelist_$enddate.log