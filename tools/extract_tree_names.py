#!/usr/bin/env python

from Bio import Phylo
from optparse import OptionParser
import sys

"""extract names from a provided tree, in order"""

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def main(tree):
    mytree = Phylo.read(tree, 'newick')
    tree_names = [ ]
    for clade in mytree.find_clades():
        if clade.name:
            print clade.name


if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage) 
    parser.add_option("-t", "--tree", dest="tree",
                      help="input tree in Newick format [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    options, args = parser.parse_args()
    
    mandatories = ["tree"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.tree)
