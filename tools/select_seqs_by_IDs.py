#!/usr/bin/python

"""takes a list of record.ids and returns to you the sequences 
from a fasta list that are part of the list"""

from Bio import SeqIO
from optparse import OptionParser
import sys

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def main(in_fasta, ids, out_fasta):
    infile = open(in_fasta, "U")
    data = open(ids, "U").read().splitlines()
    output_handle = open(out_fasta, "w")
    seqrecords=[ ]
    for record in SeqIO.parse(infile, "fasta"):
        if record.id in data:
            seqrecords.append(record)
    SeqIO.write(seqrecords, output_handle, "fasta") 
    infile.close()
    output_handle.close()

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage) 
    parser.add_option("-i", "--input_fasta", dest="in_fasta",
                    help="/path/to/input fasta [REQUIRED]",
                    action="callback", callback=test_file, type="string")
    parser.add_option("-d", "--headers", dest="ids",
                    help="/path/to/id file [REQUIRED]",
                    action="callback", callback=test_file, type="string")
    parser.add_option("-o", "--output_fasta", dest="out_fasta",
                    help="/path/to/output fasta [REQUIRED]",
                    action="store", type="string")

    options, args = parser.parse_args()
    
    mandatories = ["in_fasta", "ids", "out_fasta"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)
            
    main(options.in_fasta, options.ids, options.out_fasta)
