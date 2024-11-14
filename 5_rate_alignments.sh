#!/bin/bash

#SBATCH --partition rosa.p
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=0-10:00
#SBATCH --output=/user/rego3475/master_output/logs/5_rate_alignments.%j.out
#SBATCH --error=/user/rego3475/master_output/logs/5_rate_alignments.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ronja.roesner@uol.de

