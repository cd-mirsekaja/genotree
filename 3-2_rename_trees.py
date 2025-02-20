#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 17:10:48 2024

@author: Ronja Roesner

Script for renaming the tips of phylogenetic trees with the corresponding index from the database
"""

import argparse, sqlite3, os, re

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tree' , dest = 'tree_file' , type = str , default= None , required= True, help = 'input treefile')
parser.add_argument('-d', '--database' , dest = 'library' , type = str , default= None , required= True, help = 'input sql database')
args, unknown = parser.parse_known_args()

locus_id = args.tree_file.split("-")[0]
db_file=args.library

# remove database file if toggle is set to 1, else print different statements
if os.path.isfile(db_file):
	# establishes connection to the database and creating an empty file if there is none
	db_conn=sqlite3.connect(db_file)
	# creates new cursor object to interact with the database
	c=db_conn.cursor()
else:
	raise Exception('Database does not exist. Exiting.')


# writes contents of input tree to variable treefile_content
with open(args.tree_file, 'r') as file:
	treefile_content = file.read()

# set variables for counting and creating a list of genomes
genome_count=0
genome_list=[]

# query to get the index and accession number from the database
query="SELECT IDX, AccessionNumber FROM ids"
c.execute(query)
# iterates over the rows of the query result to replace the accession number with the index
for row in c:
	idx,acc_number=row
	replacement_str=f"{idx}"
	
	if acc_number in treefile_content:
		treefile_content = treefile_content.replace(acc_number, replacement_str)
		genome_list.append(acc_number)
		genome_count=genome_count+1

# converts the list of genomes to a string
genome_string=str(genome_list)

# replaces all unneccessary clutter in the treefile with ":"
treefile_content = re.sub(r'-[A-Z0-9_]+\.?[0-9]*-\d+-\d+:', ':', treefile_content)

# closes connection to the database
db_conn.close()

# writes treefile_content to new file
with open(locus_id + "-" + str(genome_count) + "_genomes" + "-renamed" + ".treefile","w") as tree_out:
	tree_out.write(treefile_content)
	print("Renamed treefile saved to ",tree_out.name)

with open("genome_list"+".log","a") as outfile:
	outfile.write('\n')
	outfile.write(args.tree_file+": "+str(genome_count)+" genomes")
	outfile.write('\n')
	outfile.write(genome_string)
	
	print("List of genomes saved to ",outfile.name)

