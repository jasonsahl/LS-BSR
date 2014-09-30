#!/usr/bin/env python

"""re-orders a BSR matrix, based on a user-defined
list of genomes.  This will typically be associated
with the order of genomes in a phylogeny"""

import sys
import optparse
from ls_bsr.util import transpose_matrix, reorder_matrix

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
