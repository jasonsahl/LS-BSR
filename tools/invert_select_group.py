#!/usr/bin/python

from sys import argv
import sys

try:
    with open(argv[1]) as name_file:
        names = name_file.read().splitlines()
    with open(argv[2]) as count_file:
        counts = count_file.read().splitlines()
except:
    print("usage: invert_select_group.py group1_ids.txt names.txt > group2_ids.txt")
    sys.exit()

for x in counts:
    if x not in names:
        print(x)
    else:
        pass
