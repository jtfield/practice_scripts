#! /usr/bin/python3

import os
import sys

folder = sys.argv[1]
print(folder)
for file in os.listdir(folder):
    if "-gb" in file:
        if "htm" not in file:
            print(file)
            filepath = os.path.join(folder, file)
            f = open(filepath, 'r')
            fread = f.readlines()
            for line in fread:
                if ">" in line:
                    if "PHOM" not in line:
                        if "PDEL" not in line:
                            print(line)
            f.close()
