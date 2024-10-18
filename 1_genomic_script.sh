#!/bin/bash

#SBATCH --partition rosa.p
#SBATCH --cpus-per-task=10
#SBATCH --mem=80G
#SBATCH --time=0-22:00
#SBATCH --output=./logs/automatic/1_genomic_script.%j.out
#SBATCH --error=./logs/automatic/1_genomic_script.%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ronja.roesner@uol.de

module load hpc-env/13.1
module load HMMER/3.4-gompi-2023a
module load Python/3.11.3-GCCcore-13.1.0
module load Biopython/1.81-foss-2023a

genome="$1"
cpu_count=$(expr "$2"+0)

# extract the accession number from the file name
species_id=$(echo "${genome%%.fasta}")

# terminates the script if the accession number is on the forbidden list
if grep -Fxq "$species_id" ~/master_input/forbidden_accnumbers.txt
then exit 0
fi

# display start time in log file
echo === start time for $species_id is $(date '+%H:%M:%S') ===

# make a folder for the nhmmer files of this genome
mkdir ./nhmmer-tables/$species_id-tables

for file in ./hmm/*.hmm; do
	# extract the locus id from the file name
	locus_id=$(echo "${file##*/}" | cut -d'-' -f1)
	
	echo === starting nhmmer search for $locus_id ===
	# run nhmmer for current genome
	nhmmer --cpu 8 --noali -T 100 --tblout "nhmmer-tables/$species_id-tables/$species_id-$locus_id-table.txt" "$file" "/nfs/data/zapp2497/genomes/raw_fasta/$genome"
	
	# rename nhmmer output to be better readable
	#for f in ./nhmmer-tables/$species_id-tables/$species_id-$file.table.txt; do mv -- "$f" "${f%-NoDups_OK.hmm.table.txt}-table.txt"; done

	echo === starting hit search for $locus_id ===
	# find hits for all loci in current genome
	python3 ~/Main_analysis_II/2_hmmer2fasta_script.py -i1 "nhmmer-tables/$species_id-tables/$species_id-$locus_id-table.txt" -i2 "$genome"
	mv $species_id-$locus_id-hits.fasta hits/$locus_id-allhits
done

echo === zipping and moving files ===
# zip nhmmer tables for this genome and delete empty folder (currently also zips the parent folder 'nhmmer-tables')
zip $species_id-nhmmer.zip ./nhmmer-tables/$species_id-tables/* -m
mv $species_id-nhmmer.zip ./nhmmer-tables
rm -r ./nhmmer-tables/$species_id-tables

# move all hit files to their respective subfolder
#mv $species_id-*-hits.fasta ./hits

# display end time in log file
echo === end time for $species_id is $(date '+%H:%M:%S') ===