#!/usr/bin/env python

"""If given a GenBank file,
extract the locus tags into
a new file"""

from sys import argv
from Bio import SeqIO
import sys
import os

try:
    infile = open(argv[1], "U")
except:
    print "usage: script genbank_file.gbk"
    sys.exit()

record = SeqIO.read(infile, "genbank")
name = argv[1].split(".gbk")
output_handle = open("%s.fasta" % name[0], "w")
count = 0
for feature in record.features:
    if feature.type == "gene":
        count = count + 1
        feature_name = feature.qualifiers["locus_tag"]
        feature_seq = feature.extract(record.seq)
        # Simple FASTA output without line wrapping:
        output_handle.write(">" + "".join(feature_name) + "\n" + str(feature_seq) + "\n")
output_handle.close()
print(str(count) + " CDS sequences extracted")
