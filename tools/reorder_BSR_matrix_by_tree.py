#!/usr/bin/env python

"""re-orders a BSR matrix, based on a user-defined
list of genomes.  This will typically be associated
with the order of genomes in a phylogeny"""

import sys
import optparse
from Bio import Phylo

def transpose_matrix(matrix):
    out_matrix = open("tmp.matrix", "w")
    in_matrix = open(matrix, "U")
    reduced = [ ]
    for line in in_matrix:
        fields = line.split("\t")
        reduced.append(fields)
    test=map(list, zip(*reduced))
    for x in test:
        print >> out_matrix, "\t".join(x)
    in_matrix.close()
    out_matrix.close()
    
def reorder_matrix(in_matrix, names):
    my_matrix = open(in_matrix, "U")
    outfile = open("reordered_matrix.txt", "w")
    firstLine = my_matrix.readline()
    print >> outfile, firstLine,
    my_matrix.close()
    for name in names:
         for line in open(in_matrix, "U"):
            fields = line.split("\t")
            if name == fields[0]:
                print >> outfile, line,
           
def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def parse_tree(tree):
    names = []
    mytree = Phylo.read(tree, 'newick')
    tree_names = [ ]
    for clade in mytree.find_clades():
        if clade.name:
            names.append(clade.name)
    return names

def main(matrix, tree):
    transpose_matrix(matrix)
    names = parse_tree(tree)
    reorder_matrix("tmp.matrix", names)

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
