#!/usr/bin/env python

"""calculate several pan-genome type stats
from a BSR matrix"""

from optparse import OptionParser
from ls_bsr.util import get_core_gene_stats
from ls_bsr.util import get_frequencies
import sys

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def main(matrix, threshold, lower):
    get_core_gene_stats(matrix, threshold, lower)
    get_frequencies(matrix, threshold)
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage) 
    parser.add_option("-b", "--bsr_matrix", dest="matrix",
                      help="/path/to/bsr_matrix [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-t", "--threshold", dest="threshold",
                      help="upper threshold for ORF presence, defaults to 0.8",
                      action="store", default="0.8", type="float")
    parser.add_option("-l", "--lower", dest="lower",
                      help="lower threshold for ORF presence, defaults to 0.4",
                      action="store", default="0.4", type="float")
    options, args = parser.parse_args()
    
    mandatories = ["matrix"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.matrix, options.threshold, options.lower)
