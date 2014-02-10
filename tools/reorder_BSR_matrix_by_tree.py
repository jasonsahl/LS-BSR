#!/usr/bin/env python

"""re-orders a BSR matrix, based on a user-defined
list of genomes.  This will typically be associated
with the order of genomes in a phylogeny"""

import os
import sys
import optparse

from ls_bsr.util import transpose_matrix
from ls_bsr.util import reorder_matrix
from ls_bsr.util import parse_tree

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def main(matrix, tree):
    transpose_matrix(matrix)
    names = parse_tree(tree)
    reorder_matrix("tmp.matrix", names)
    os.system("rm tmp.matrix")

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-b", "--matrix", dest="matrix",
                      help="/path/to/BSR matrix [REQUIRED]",
                      type="string", action="callback", callback=test_file)
    parser.add_option("-t", "--tree", dest="tree",
                      help="/path/to/tree in Newick format [REQUIRED]",
                      type="string", action="callback", callback=test_file)

    options, args = parser.parse_args()

    mandatories = ["matrix","tree"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.matrix, options.tree)
