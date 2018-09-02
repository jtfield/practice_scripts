#! /usr/bin/python3

import sys
import re

infile1 = sys.argv[1]
#infile2 = sys.argv[2]
#outfile = sys.argv[3]

#output = open("uncontam_dupe_hcolv5.fasta",'w')

trans_file = open(infile1).read().split(">")
clean_trans = open("clean_hcolv5.fasta","w")

regexes = ['hcolumbae_DN110704_c0_g1_i1',
'hcolumbae_DN110704_c0_g2_i1',
'hcolumbae_DN46231_c0_g1_i1',
'hcolumbae_DN46231_c0_g2_i1',
'hcolumbae_DN61095_c0_g2_i1',
'hcolumbae_DN62752_c0_g1_i1',
'hcolumbae_DN66885_c0_g1_i1',
'hcolumbae_DN66885_c0_g1_i2',
'hcolumbae_DN71510_c0_g2_i1',
'hcolumbae_DN76803_c2_g2_i2',
'hcolumbae_DN80074_c1_g1_i1',
'hcolumbae_DN80074_c1_g3_i1',
'hcolumbae_DN80074_c1_g3_i2',
'hcolumbae_DN80861_c0_g1_i1',
'hcolumbae_DN80861_c0_g2_i1',
'hcolumbae_DN92241_c0_g1_i1',
'hcolumbae_DN92241_c0_g2_i1']

combined = "(" + ")|(".join(regexes) + ")"

for chunk in trans_file:
    if re.match(combined, chunk):
        print("Some regex matched!")
    else:
        print("nothing matched")
        clean_trans.write(">")
        clean_trans.write(chunk)
