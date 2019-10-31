#! /usr/bin/python3

import sys
import re

# USESAGE: this script will take in a .fasta file, a file containing a list of sequence names you wish to remove
# and the name of the output .fasta file you wish to produce
# example command: python3 /path/to/contam_remove.py /path/to/fasta_file /path/to/sequence_to_remove_file /path/to/output_fasta_file

# initialize user supplied arguments
infile1 = sys.argv[1]
infile2 = sys.argv[2]
outfile = sys.argv[3]

# open FASTA file for reading
trans_file = open(infile1).read().split(">")

# open output file for writing
clean_trans = open(outfile,"w")

# open and parse contamination file
contam_list = []
contam_file = open(infile2).read().split()
for line in contam_file:
    contam_list.append(line)

print("Lists of contaminated sequences created")
print("searching and removing contamination sequences from contig file")

combined = "(" + ")|(".join(contam_list) + ")"

for chunk in trans_file:
    if re.match(combined, chunk):
        print("Some regex matched!")
    else:
        #print("nothing matched")
        clean_trans.write(">")
        clean_trans.write(chunk)

clean_trans.close()
