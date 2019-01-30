#!/usr/bin/env python

"""If given a GenBank file,
extract the locus tags into
a new file"""

from sys import argv
from Bio import SeqIO
import sys
import os

try:
    name = argv[1].split(".gbk")
    output_handle = open("%s.fasta" % name[0], "w")
    count = 0
    with open(argv[1]) as infile:
        for record in SeqIO.parse(infile, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    try:
                        feature_name = feature.qualifiers["locus_tag"]
                        feature_seq = feature.extract(record.seq)
                        feature_product = feature.qualifiers["product"]
                        count = count + 1
                        # Simple FASTA output without line wrapping:
                        output_handle.write(">" + "".join(feature_name) + "|" + "".join(feature_product) + "\n" + str(feature_seq) + "\n")
                    except:
                        pass
    output_handle.close()
except:
    print("usage: script genbank_file.gbk")
    sys.exit()
if count == 0:
    print("Problem with GenBank file. Make sure that it contains locus_tags")
else:
    print(str(count) + " CDS sequences extracted")
