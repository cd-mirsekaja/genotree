#!/bin/bash

#SBATCH --partition rosa.p
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=0-48:00
#SBATCH --output=/user/rego3475/master_output/logs/5_process_alignments.%j.out
#SBATCH --error=/user/rego3475/master_output/logs/5_process_alignments.%j.err
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

# use for testing
alignment_count=2

# use for full dataset
#alignment_count=$(ls ~/master_input/all_hits_aligned/ | wc -l)


alignment_files=$(ls ~/master_input/all_hits_aligned/ | head -n $alignment_count)

echo === starting alignment scoring at $(date '+%d.%m.%Y %H:%M:%S') ===
# prepares a function that runs the script for scoring each alignment
scoring_function() {
	locus_in=$(echo "${1%%-*}")
	sbatch ~/genotree/6_rate_alignments.sh $locus_in
}

# exports the function so GNU Parallel can access it
export -f scoring_function

# runs 6_rate_alignments.sh in parallel for each specified exon
echo "$alignment_files" | parallel --eta --jobs $alignment_count scoring_function

# waits for all genome files to be processed and nhmmer-tables to be created by checking the squeue command for a certain term
while squeue -u $USER | grep -q "6_"; do wait; done

