#! /usr/bin/python3

import sys

infile = sys.argv[1]
nucleotide_num = sys.argv[2]

num = int(nucleotide_num)

data = open(infile).read().split(">")

for line in data:
    print(">" + line[0:num])
