#!/bin/bash

#SBATCH --partition rosa.p
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=0-52:00
#SBATCH --output=/user/rego3475/master_output/logs/0_job_main_script.%j.out
#SBATCH --error=/user/rego3475/master_output/logs/0_job_main_script.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ronja.roesner@uol.de,simon.kaefer@uol.de

# get the starting time
startdate=$(date '+%Y_%m_%d-%H_%M_%S')
echo === start date and time is $startdate ===

# load necessary modules
module load hpc-env/13.1
module load parallel/20230822-GCCcore-13.1.0
module load HMMER/3.4-gompi-2023a
module load Python/3.11.3-GCCcore-13.1.0
#module load R/4.3.1-foss-2023a

# make working directory and move into it
mkdir $WORK/wd-$startdate
cd $WORK/wd-$startdate

# make output directories
mkdir hits hits-aligned hits-combined hmm nhmmer-tables treefiles-original treefiles-renamed treefiles-final trash
mkdir logs
mkdir logs/automatic
mkdir logs/manual

# make log files
touch logfile.log
touch genome_list.log

# when testing, use these - sets the amount of genomes and loci for checking
#genomecount=10
#locuscount=10

# when running with full dataset, use these
genomecount=$(ls /nfs/data/zapp2497/genomes/raw_fasta/ | wc -l)
locuscount=$(ls ~/master_input/locus-fa_1105_nucleotide/ | wc -l)

# sets the amount of cpus that nhmmer should use
cpu_count=8

# add the working directory and start time to the logfile
echo Working Directory: $(pwd) >> logfile.log
(echo  ;echo Starting at: $startdate) >> logfile.log

# declares a variable genomic_files containing the names all genomes in the dataset
genomic_files=$(ls /nfs/data/zapp2497/genomes/raw_fasta/| sort -R | head -n $genomecount)

# declares a variable locus_files containing the names of all exons in the dataset
locus_files=$(ls ~/master_input/locus-fa_1105_nucleotide/ | head -n $locuscount)

echo === making hmm files at $(date '+%H:%M:%S') ===

# build .hmm files for all .al files in a subfolder of the working directory
for file in $locus_files; do 
	locus_id=$(echo "${file%%-*}")
	mkdir ./hits/$locus_id-allhits
	hmmbuild $file.hmm ~/master_input/locus-fa_1105_nucleotide/$file
done

#remove unnecessary clutter in file file names by renaming all .hmm files
for f in *.fa.hmm; do mv -- "$f" "${f%.fa.hmm}.hmm"; done

# moves all .hmm files to their respective subdirectory
mv *.hmm ./hmm

echo === starting genomic analysis scripts at $(date '+%d.%m.%Y %H:%M:%S') ===
# prepares a function that runs the script for creating hit.fasta files for every genome
genomic_function() {
	genomic_in="$1"
	sbatch ~/Main_analysis_II/1_genomic_script.sh $genomic_in $cpu_count
}

# exports the function so GNU Parallel can access it
export -f genomic_function

# runs genomic_script.sh in parallel for each specified genome
echo "$genomic_files" | parallel --eta --jobs $genomecount genomic_function

# waits for all genome files to be processed and nhmmer-tables to be created by checking the squeue command for a certain term
while squeue -u $USER | grep -q "1_"; do wait; done

echo === starting locus analysis scripts at $(date '+%d.%m.%Y %H:%M:%S') ===
# prepares a function that runs the script for creating trees for each locus
align_tree_function() {
	locus_in=$(echo "${1%%-*}")
	sbatch ~/Main_analysis_II/3_align_tree_script.sh $locus_in
}

# exports the function so GNU Parallel can access it
export -f align_tree_function

# runs align_tree_script.sh in parallel for each specified locus
echo "$locus_files" | parallel --eta --jobs $locuscount align_tree_function

# waits for all trees to be created by checking the squeue command for a certain term
while squeue -u $USER | grep -q "3_"; do wait; done

# add a list of all used genomes to the logfile including the number of times a genome was found in the trees
(echo  ;echo $genomecount genomes: ) >> logfile.log
tree_count=$(ls treefiles-renamed/*.treefile | wc -l)
for file in $genomic_files;do
	acc_number=$(echo "${file%%.fasta}")
	acc_count=$(grep $acc_number genome_list.log | wc -l)	
	(echo -n "$acc_number | ";echo -n $(python3 ~/Main_analysis_II/5_speciesinfo.py -n "${file%%.fasta}" -x ~/master_input/genome_master_library.xlsx);echo " | found in $acc_count out of $tree_count trees") >> logfile.log

done


# add a list of all loci to the logfile
(echo  ;echo $locuscount loci: )>> logfile.log
for file in $locus_files; do
	locus_id=$(echo "${file%%-NoDups_OK.fa}")
	locus_line=$(grep $locus_id genome_list.log)
	loc_sp_count=$(echo "${locus_line#*: }")
	(echo "$locus_id | $loc_sp_count in tree") >> logfile.log
done

# make a combined tree out of individual gene trees and run astral on it (program installed locally)
cat treefiles-renamed/*.treefile > treefiles-final/all-loci_combined.treefile
~/programs/ASTER-Linux/bin/astral-pro3 -t 8 -o treefiles-final/all-loci_astralpro3.treefile treefiles-final/all-loci_combined.treefile 2>astralpro3.log
~/programs/ASTER-Linux/bin/astral4 -t 8 -o treefiles-final/all-loci_astral4.treefile treefiles-final/all-loci_combined.treefile 2>astral4.log

enddate=$(date '+%Y_%m_%d-%H_%M_%S')

echo === moving files ===
# make folder for all outputs of this run
mkdir ~/master_output/$enddate-all_out_FULL-DATASET

# add the ending time to the logfile and move all log files
(echo  ;echo Ending at: $enddate) >> logfile.log
mv logfile.log logs/manual/logfile_$enddate.log
mv genome_list.log logs/manual/genomelist_$enddate.log
mv astralpro3.log logs/automatic/astral3pro_$enddate.log
mv astral4.log logs/automatic/astral4_$enddate.log

# move all output folders
mv * ~/master_output/$enddate-all_out_FULL-DATASET

cd ~/Main_analysis_II/
rm -r $WORK/wd-$startdate

echo === end date and time is $enddate ===