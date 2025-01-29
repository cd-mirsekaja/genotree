#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 17:39:48 2025

@author: Ronja Rösner

Gets a aligroove score matrix, an alignment file, a maximum threshold and a locus id as inputs. Outputs table containing all values from the input matrix that are below the threshold with their corresponding row and column names as well as an alignment file filtered for the threshold.
"""

import pandas as pd
import argparse


# add input arguments
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--input_table' , dest='input_table' , type=str , default=None , required=True, help='input aligroove score matrix')
parser.add_argument('-a', '--input_alignment' , dest='input_alignment' , type=str , default=None , required=True, help='input renamed alignment file')
parser.add_argument('-l', '--locus_id', dest='locus_id' , type=str, default=None , required=True, help='input locus id')
parser.add_argument('-t', '--threshold', dest='threshold', type=float, default=0.2, required=False, help='threshold value (default: 0.2)')
args, unknown = parser.parse_known_args()


# import matrix with the total mean and median scores for all loci
#total_score_matrix = pd.read_csv(f"{args.dir_scores}/total_scores.csv",sep=";")

# import matrix with all values for this locus from file
score_matrix = pd.read_csv(args.input_table,sep="\t",index_col=0)
# create filtered matrix with only values lower than or equal to the threshold
filtered_matrix = score_matrix.where(score_matrix <= args.threshold)

# the following 4 lines were written with the help of ChatGPT
# Stack the filtered matrix to get row indices, column names, and values
low_values = filtered_matrix.stack().reset_index()
low_values.columns = ['Row', 'Column', 'Value']
low_values['Value'] = low_values['Value'].round(2)  # Round values to 2 decimals

# exit the program if no values are below the threshold
if low_values.empty:
	print(f"No values found below the threshold ({args.threshold}). Exiting.")
	exit(0)

# save low_values to file
with open(f"{args.locus_id}_{args.threshold}-values"+".csv","w") as file:
	low_values.to_csv(file, sep=';', index=False, header=['Column', 'Row', 'Value (rounded)'])
	print(f"Low values saved as {file.name}")


alignment_content=open(args.input_alignment)
new_alignment=open(f"{args.locus_id}_{args.threshold}-filtered.fasta","w")

# the following 12 lines were written with the help of ChatGPT
write_flag = True  # Determines whether to write lines to output
for line in alignment_content:
	write_flag = True  # Determines whether to write lines to output
	for line in alignment_content:
		if line.startswith(">"):  # If it's a headline
			item_id = line[1:].strip()  # Extract item ID (remove > and strip whitespace)
			if item_id in low_values:  # If headline matches an item in the table
				write_flag = False  # Stop writing until next headline
			else:
				write_flag = True  # Resume writing for other headlines
		if write_flag:
			new_alignment.write(line)  # Write line to the output file

print(f"Filtered alignment written to {new_alignment.name}")

alignment_content.close()
new_alignment.close()


