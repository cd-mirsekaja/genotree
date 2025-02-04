#!/bin/bash

#SBATCH --partition rosa.p
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=0-50:00
#SBATCH --output=/user/rego3475/master_output/logs/2_process_alignments.%j.out
#SBATCH --error=/user/rego3475/master_output/logs/2_process_alignments.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ronja.roesner@uol.de

# get the starting time
startdate=$(date '+%Y_%m_%d-%H_%M_%S')
echo === start date and time is $startdate ===

# load necessary modules
module load hpc-env/13.1
module load parallel/20230822-GCCcore-13.1.0
module load Python/3.11.3-GCCcore-13.1.0

mkdir $WORK/wd-al_scoring-$startdate
cd $WORK/wd-al_scoring-$startdate

mkdir logs
mkdir aligroove_output alignments_filtered filter_tables
mkdir aligroove_output/txt aligroove_output/svg 

# low amount of files, use for testing
#alignment_count=5

# all files, use for full dataset
alignment_count=$(ls ~/master_input/all_hits_aligned_renamed/*-renamed.fasta | wc -l)

# save all relevant alignment files in a variable
alignment_files=$(ls ~/master_input/all_hits_aligned_renamed/*-renamed.fasta | head -n $alignment_count)

echo === starting alignment scoring at $(date '+%d.%m.%Y %H:%M:%S') ===
# prepares a function that runs the script for scoring each alignment
scoring_function() {
	locus_in=$(echo "${1##*/}" | cut -d'-' -f1)
	sbatch ~/genotree/2-1_rate_alignments.sh $locus_in
}

# exports the function so GNU Parallel can access it
export -f scoring_function

# runs 2-1_rate_alignments.sh in parallel for each specified exon
echo "$alignment_files" | parallel --eta --jobs $alignment_count scoring_function

# waits for all alignments to be rated by checking the squeue command for a certain term
while squeue -u $USER | grep -q "2-1"; do sleep 1; wait; done

echo === finished alignment scoring at $(date '+%d.%m.%Y %H:%M:%S') ===

mv ~/master_input/all_hits_aligned_renamed/*-renamed.txt ./aligroove_output/txt
mv ~/master_input/all_hits_aligned_renamed/*-renamed.svg ./aligroove_output/svg

echo === calculating total mean and median scores ===
printf "locusID;scoreMedian;scoreMean\n" > total_scores.csv

# calculate total mean and median scores for all alignments
for file in ./aligroove_output/txt/*-renamed.txt; do
    locus_id=$(echo "${file##*/}" | cut -d'_' -f5 | cut -d'-' -f1)
    echo - $locus_id -
    python3 ~/genotree/utilities/total_score.py -f $file -l $locus_id
done

# move score table to output folder
mv total_scores.csv ./aligroove_output

# set filter threshold
threshold=0.35
echo === filtering alignments for score threshold $threshold ===
for file in $alignment_files; do
    # get ID for current locus
    locus_id=$(echo "${file##*/}" | cut -d'-' -f1)
    aligroove_matrix=./aligroove_output/txt/AliGROOVE_seqsim_matrix_$locus_id-renamed.txt
    echo Locus ID: $locus_id
    # filters alignments for current locus for all exons with a score under the given threshold
    python3 ~/genotree/2-2_filter_alignments.py -d $aligroove_matrix -a $file -l $locus_id -t $threshold
done

echo === moving files ===

mv *-filtered.fasta ./alignments_filtered
mv *-values.csv ./filter_tables


enddate=$(date '+%Y_%m_%d-%H_%M_%S')
mkdir ~/master_output/refined_alignments/$enddate-alignments_FULL-DATASET
mv * ~/master_output/refined_alignments/$enddate-alignments_FULL-DATASET

cd ~/genotree
#rm -r $WORK/wd-al_scoring-$startdate