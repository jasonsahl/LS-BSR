#!/usr/bin/env python

"""If given a GenBank file,
extract the locus tags into
a new file"""

from __future__ import print_function
from sys import argv
from Bio import SeqIO
import sys
import os

try:
    infile = open(argv[1], "U")
except:
    print("usage: script genbank_file.gbk")
    sys.exit()

curr_dir = os.getcwd()
curr_path=os.path.abspath("%s" % curr_dir)
record = SeqIO.read(infile, "genbank")
name = argv[1].split(".gbk")
output_handle = open("%s/%s.fasta" % (curr_path,name[0]), "w")
count = 0
for feature in record.features:
    if feature.type == "gene":
        count = count + 1
        try:
            feature_name = feature.qualifiers["locus_tag"]
            feature_seq = feature.extract(record.seq)
            print(str(record.id)+"\t"+"".join(feature_name)+"\t"+str(feature.location))
            # Simple FASTA output without line wrapping:
            output_handle.write(">"+"".join(feature_name)+"\n"+str(feature_seq)+"\n")
        except:
            print("problem encountered extracting: %s, skipping.." % "".join(feature_name))
output_handle.close()
os.system("mv %s.fasta %s" % (name[0],curr_path))
print(str(count) + " CDS sequences extracted")
