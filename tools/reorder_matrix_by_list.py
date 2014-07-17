#!/usr/bin/env python

"""re-orders a BSR matrix, based on a user-defined
list of genomes.  This will typically be associated
with the order of genomes in a phylogeny"""

import sys
import optparse

def transpose_matrix(matrix):
    out_matrix = open("tmp.matrix", "w")
    in_matrix = open(matrix, "U")
    reduced = [ ]
    for line in in_matrix:
        newline = line.strip("\n")    
        fields = newline.split("\t")
        reduced.append(fields)
    test=map(list, zip(*reduced))
    for x in test:
        print >> out_matrix, "\t".join(x)
    in_matrix.close()
    out_matrix.close()
    
def reorder_matrix(in_matrix, genomes):
    my_matrix = open(in_matrix, "U")
    outfile = open("reordered_matrix.txt", "w")
    firstLine = my_matrix.readline()
    print >> outfile, firstLine,
    for stuff in open(genomes, "U").read().splitlines():
        for line in open(in_matrix, "U"):
	    newline = line.strip()        
            fields = newline.split("\t")
            if stuff == fields[0]:
                print >> outfile, line,
           
def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def main(matrix, genomes):
    transpose_matrix(matrix)
    reorder_matrix("tmp.matrix", genomes)

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-b", "--matrix", dest="matrix",
                      help="/path/to/BSR matrix [REQUIRED]",
                      type="string", action="callback", callback=test_file)
    parser.add_option("-g", "--genomes", dest="genomes",
                      help="/path/to/new_line delimited genomes file [REQUIRED]",
                      type="string", action="callback", callback=test_file)

    options, args = parser.parse_args()

    mandatories = ["matrix","genomes"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.matrix, options.genomes)
