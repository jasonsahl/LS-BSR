#!/usr/bin/env python

"""from a BSR matrix, filter
out un-desired genomes"""

from optparse import OptionParser
from collections import deque
import sys
from ls_bsr.util import filter_genomes
from ls_bsr.util import filter_matrix

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def main(in_matrix, prefix, genomes):
    to_keep=filter_genomes(genomes, in_matrix)
    filter_matrix(to_keep, in_matrix, prefix)
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage) 
    parser.add_option("-b", "--in_matrix", dest="in_matrix",
                      help="/path/to/BSR matrix [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-p", "--out_prefix", dest="prefix",
                      help="output naming prefix [REQUIRED]",
                      action="store", type="string")
    parser.add_option("-g", "--genomes", dest="genomes",
                      help="/path/to/genomes_file (to remove) (new line delimited) [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    options, args = parser.parse_args()
    
    mandatories = ["in_matrix", "prefix", "genomes"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.in_matrix,options.prefix,options.genomes)
