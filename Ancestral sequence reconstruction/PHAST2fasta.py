#!/cluster/apps/gdc/python/3.6.1/bin/python3.6
### This script converts the output of PHAST prequel to fasta format ###
## It does so by calling bases according the maximum probability that the based was the ancestral base (from prequel) and by applying a probability threshold below which the base is returned as N.

# Import modules
import os
import pandas
import argparse

# Change working directory
#os.chdir("/Users/luqman/Desktop/Ancestral sequence reconstruction")

## We parse arguments from the command line using argparse 
parser = argparse.ArgumentParser(description="PHAST prequel to fasta (Luqman 2018)", add_help=True)
parser.add_argument("filename", action="store", help="The input file (PHAST prequel output)")
parser.add_argument("-p", action="store", default="header", dest="fasta_header", help="Fasta header")
parser.add_argument("-o", action="store", default="scaffold", dest="output", help="Prefix for output files.")
args = parser.parse_args()

filename = args.filename
fasta_header = args.fasta_header
output_prefix = args.output

# Import data
scaffold_probs_raw = pandas.read_table(filename)
#scaffold_probs_raw = pandas.read_table("scaffold_2_anc.CAGATC-Longicaulis.probs")
scaffold_probs_raw.columns = ['A', 'C', 'G', 'T']

# Make copy
scaffold_probs = scaffold_probs_raw.copy()

# If maximum probability < 0.75, output N, else X
def max_threshold(A,C,G,T):
	max_ACTG = max([A,C,G,T])
	if max_ACTG < 0.75:
		base = "N"
	else:
		base = "X"
	return base

scaffold_probs["max_pass"] = scaffold_probs.apply(lambda x: max_threshold(x["A"], x["C"], x["G"], x["T"]),axis=1)

# Find base with max probability
scaffold_probs["max"] = scaffold_probs_raw[["A","C", "G", "T"]].idxmax(axis=1)

# Output base with max probability given that it passes a given threshold
def base_call(max_pass,max):
	if max_pass == "X":
		base = max
	elif max_pass == "N":
		base = max_pass
	return base

scaffold_probs["base"] = scaffold_probs.apply(lambda x: base_call(x["max_pass"], x["max"]),axis=1)

# Output into fasta format
# Convert series to list
called_bases_list = scaffold_probs["base"].tolist()
# Convert list to string
called_bases_string = ''.join(called_bases_list)
# Split string into 60 character chunks
n=60
called_bases_fasta_list = [called_bases_string[i:i+n] for i in range(0, len(called_bases_string), n)]


# Write to file
with open(output_prefix + ".fasta", 'w') as fasta_file:
	fasta_file.writelines(">" + fasta_header + "\n")
	#fasta_file.writelines(">" + "header_1" + "\n")
	fasta_file.writelines("%s\n" % l for l in called_bases_fasta_list)
	
