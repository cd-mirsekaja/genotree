#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 17:10:48 2024

@author: Ronja Roesner
"""

import re
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tree' , dest = 'tree_file' , type = str , default= None , required= True, help = 'input treefile')
parser.add_argument('-x', '--table' , dest = 'overview_library' , type = str , default= None , required= True, help = 'input excel table')
args, unknown = parser.parse_known_args()

# imports reference table as dataframe table
table =  pd.read_excel(args.overview_library, header=0)
locus_id = args.tree_file.split("-")[0]

# writes contents of input tree to variable treefile_content
with open(args.tree_file, 'r') as file:
	treefile_content = file.read()

genome_count=0
genome_list=[]

# iterates through all rows in the table and replaces the accession numbers with the replacement_string
for index, row in table.iterrows():
	spec_id=row['Index']
	accession_number=row['AccessionNumber']
	
	replacement_string = f"{spec_id}"
	
	if accession_number in treefile_content:
		treefile_content = treefile_content.replace(accession_number, replacement_string)
		genome_list.append(accession_number)
		genome_count=genome_count+1
	

species_string=str(genome_list)

# replaces all unneccessary clutter with ":"
treefile_content = re.sub(r'\|locus_E\d+-NoDups_OK\|[^:]+:', ':', treefile_content)


# writes treefile_content to new file
tree_out = open(locus_id + "-" + str(genome_count) + "_genomes" + "-renamed" + ".treefile","w")
tree_out.write(treefile_content)
tree_out.close()

genome_list_out = open("genome_list"+".log","a")
genome_list_out.write('\n')
genome_list_out.write(args.tree_file+": "+str(genome_count)+" genomes")
genome_list_out.write('\n')
genome_list_out.write(species_string)
genome_list_out.close()

