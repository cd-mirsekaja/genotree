#!/bin/bash

#SBATCH --partition rosa.p
#SBATCH --cpus-per-task=6
#SBATCH --mem=60G
#SBATCH --time=0-12:00
#SBATCH --output=./logs/automatic/3_align_tree_script.%j.out
#SBATCH --error=./logs/automatic/3_align_tree_script.%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ronja.roesner@uol.de

module load hpc-env/13.1
module load MAFFT/7.505-GCC-13.1.0-with-extensions
module load IQ-TREE/2.2.2.7-gompi-2023a
module load Python/3.11.3-GCCcore-13.1.0

locus_id="$1"

# display start time in log file
echo === start time for $locus_id is $(date '+%H:%M:%S') ===

# combines all created .fasta files for this locus into one file
cat hits/$locus_id-allhits/* >> $locus_id-combined_hits.fasta

# zip the original hit files for this locus and remove the empty folder
zip $locus_id-allhits.zip hits/$locus_id-allhits/* -m
rm -r hits/$locus_id-allhits
mv $locus_id-allhits.zip ./hits

echo === aligning hits for $locus_id ===
# runs mafft on all combined hits for this exon
mafft --thread 6 --auto $locus_id-combined_hits.fasta > $locus_id-aligned.fasta

echo === generating tree for $locus_id ===
# creates a tree for this aligned hitfile
iqtree2 -s $locus_id-aligned.fasta -T 6 --tbe --alrt 10000

# renames tree branches to simplify analysis
python3 ~/genotree/4_treenaming_script.py -t $locus_id-aligned.fasta.treefile -x ~/master_input/genome_master_library.xlsx

echo === moving files ===
# move all combined, tree and miscellaneous files to their respective subfolders
mv $locus_id-aligned.fasta ./hits-aligned
mv $locus_id-combined_hits.fasta ./hits-combined
mv $locus_id-aligned.fasta.treefile ./treefiles-original
mv $locus_id*-renamed.treefile ./treefiles-renamed
mv $locus_id*.fasta.* ./trash

# display end time in log file
echo === end time for $locus_id is $(date '+%H:%M:%S') ===