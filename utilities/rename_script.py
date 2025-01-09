import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file' , dest = 'file' , type = str , default= None , required= True, help = 'input treefile')
args, unknown = parser.parse_known_args()

locus_id = args.file.split("-")[0]
ending = args.file.split(".")[-1]


with open(args.file, 'r') as file:
	file_content = file.read()

if ending=="fasta":
    file_content = re.sub(r'\|locus_E\d+-NoDups_OK\|[^\n]+\n', '\n', file_content)
elif ending=="treefile":
    file_content = re.sub(r'\|locus_E\d+-NoDups_OK\|[^:]+:', ':', file_content)
else:
    raise TypeError("file type not recognized")

file_out = open(locus_id+"-renamed." + ending,"w")
file_out.write(file_content)
file_out.close()