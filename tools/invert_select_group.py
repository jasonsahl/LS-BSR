#!/usr/bin/python

#usage script group1.txt names.txt > group2.txt

from sys import argv
import sys

try:
    name_file = open(argv[1], "rU")
    count_file = open(argv[2], "rU")
except:
    print "usage: invert_select_group.py group1_ids.txt names.txt > group2_ids.txt"
    sys.exit()

counts = count_file.read().splitlines()
names = name_file.read().splitlines()

for x in counts:
    if x not in names:
        print x
    else:
        pass
    
name_file.close()
count_file.close()
