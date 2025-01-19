#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 16:13:16 2024

@author: Ronja Roesner
"""

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-n', '--acnumber' , dest = 'ac_number' , type = str , default= None , required= True, help = 'input accession number')
parser.add_argument('-x', '--table' , dest = 'library' , type = str , default= None , required= True, help = 'input excel table')
args, unknown = parser.parse_known_args()

table_taxonomy =  pd.read_excel(args.library, usecols="A:E",header=0)
table_habitats = pd.read_excel(args.library, usecols="B,G:J",header=0)

def find_species(table,ac_number):
	# Search for the accession number in the table
	matched_lines = table[table[table.columns[1]] == ac_number].index.tolist()

	# output the species and taxon group
	index_values = table.iloc[matched_lines,0].values.tolist()
	species_values = table.iloc[matched_lines,2].values.tolist()
	authority_values = table.iloc[matched_lines,3].values.tolist()
	taxgroup_values = table.iloc[matched_lines,4].values.tolist()
	
	# check if anything was found in the table
	if species_values==[]:
		return "not found", "not found", "not found", "not found"
	else:
		# extract the first list element and convert it into a string
		species_str=''.join(map(str, species_values[0]))
		taxgroup_str = ''.join(map(str, taxgroup_values[0]))
		authority_str = ''.join(map(str, authority_values[0]))
		index_str = str(index_values[0])
		return species_str, authority_str, taxgroup_str, index_str

def find_habitat(table,ac_number):
	matched_lines = table[table[table.columns[0]] == ac_number].index.tolist()
	out_list=[]
	
	if table.iloc[matched_lines,1].values==1:
		out_list.append("marine")
		
	if table.iloc[matched_lines,2].values==1:
		out_list.append("brackish")
		
	if table.iloc[matched_lines,3].values==1:
		out_list.append("freshwater")
		
	if table.iloc[matched_lines,4].values==1:
		out_list.append("terrestrial")
		
	if out_list!=[]:
		out_string=', '.join(map(str,out_list))
	else:
		out_string=''.join("habitats unknown")
		
	return out_string

species,authority,taxgroup,index=find_species(table_taxonomy,args.ac_number)
habitats=find_habitat(table_habitats, args.ac_number)

#check if species was found in the reference table, print it out if yes
if species=="not found":
	out_str="species not in table"
else:
	out_str=("Index "+index+" | "+species+" "+authority+" - "+taxgroup+" ("+habitats+")")

print(out_str)
