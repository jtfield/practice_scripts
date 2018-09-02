#! /usr/bin/python3

import sys
import re

infile1 = sys.argv[1]
infile2 = sys.argv[2]
outfile = sys.argv[3]

#output = open("uncontam_dupe_hcolv5.fasta",'w')

trans_file = open(infile1).read().split(">")
clean_trans = open(outfile,"w")

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
