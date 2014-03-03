#!/usr/env/bin python

"""from a BSR matrix, report back
the number of unique CDS regions
for each genome, sorted by a given
tree"""

from optparse import OptionParser
import sys, os
from Bio import Phylo

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def get_uniques(matrix, threshold):
    in_matrix = open(matrix, "U")
    outfile = open("summary_stats.tmp.txt", "w")
    firstLine = in_matrix.readline()
    firstFields = firstLine.split()
    my_dict={}
    for line in in_matrix:
        hits=[]
        fields=line.split()
        for x in fields[1:]:
            if float(x)>=float(threshold):
                hits.append(fields.index(x))
            if float(x)>=0.40:
                hits.append("1")
        if len(hits)==2:
            for x in fields[1:]:
                if float(x)>=float(threshold):
                    try:
                        my_dict[firstFields[fields.index(x)-1]].append(fields[0])
                    except:
                        my_dict[firstFields[fields.index(x)-1]] = [fields[0]]
    for k,v in my_dict.iteritems():
        print >> outfile, k, len(v)
    for firstField in firstFields:
        if firstField not in my_dict:
            print >> outfile, firstField, "0"
        
def sort_uniques_by_tree(summary, tree):
    outfile = open("uniques_sorted_by_tree.txt", "w")
    mytree = Phylo.read(tree, 'newick')
    tree_names = [ ]
    for clade in mytree.find_clades():
        if clade.name:
            tree_names.append(clade.name)
    for tree_name in tree_names:
        for line in open(summary, "U"):
            fields = line.split()
            if fields[0] == tree_name:
                print >> outfile, line,
            else:
                pass
    
def main(matrix, tree, threshold):
    get_uniques(matrix, threshold)
    sort_uniques_by_tree("summary_stats.tmp.txt", tree)
    os.system("rm summary_stats.tmp.txt")
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-b", "--bsr_matrix", dest="matrix",
                      help="/path/to/bsr_matrix [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-r", "--tree", dest="tree",
                      help="/path/to/tree [REQUIRED]",
                      action="store", type="string")
    parser.add_option("-t", "--threshold", dest="threshold",
                      help="lower threshold for ORF presence, defaults to 0.8",
                      action="store", default="0.8", type="float")
   
    options, args = parser.parse_args()
    
    mandatories = ["matrix","tree"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.matrix, options.tree, options.threshold)
