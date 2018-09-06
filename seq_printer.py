#! /usr/bin/python3

import sys

infile1 = sys.argv[1]
number_of_seqs_to_print = sys.argv[2]

file = open(infile1).read().split(">")

counter = 0
for chunk in file:
    if counter < number_of_seqs_to_print:
        print(">" + chunk)
        counter+=1
