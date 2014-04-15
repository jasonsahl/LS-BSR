#!/usr/bin/env python

"""compares BSR values between two groups in a BSR matrix
Numpy and BioPython need to be installed.  Python version must be at
least 2.7 to use collections"""

from optparse import OptionParser
import subprocess
from ls_bsr.util import prune_matrix
from ls_bsr.util import compare_values
from ls_bsr.util import find_uniques
import sys
import os

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def add_headers(infile, outfile, lower, upper):
    file_out = open(outfile, "w")
    print >> file_out,"marker"+"\t"+"group1_mean"+"\t"+">="+str(upper)+"\t"+"total_in_group_1"+"\t"+">="+str(lower)+"\t"+"group2_mean"+"\t"+">="+str(upper)+"\t"+"total_in_group2"+"\t"+">="+str(lower)
    for line in open(infile, "U"):
        print >> file_out, line,
    file_out.close()

def main(matrix,group1,group2,fasta,upper,lower):
    prune_matrix(matrix,group1,group2)
    compare_values("group1_pruned.txt","group2_pruned.txt",upper,lower)
    subprocess.check_call("paste group1_out.txt group2_out.txt > groups_combined.txt", shell=True)
    find_uniques("groups_combined.txt",fasta)
    add_headers("groups_combined.txt","groups_combined_header.txt",lower,upper)
    os.system("rm group1_out.txt group2_out.txt")
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage) 
    parser.add_option("-b", "--bsr_matrix", dest="matrix",
                      help="/path/to/bsr_matrix [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-f", "--fasta", dest="fasta",
                      help="/path/to/ORF_fasta_file [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-1", "--group_1_ids", dest="group1",
                      help="new line separated file with group1 ids [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-2", "--group_2_ids", dest="group2",
                      help="new line separated file with group2 ids [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-u", "--upper_bound", dest="upper",
                      help="upper bound for BSR comparisons, defaults to 0.8",
                      default="0.8", type="float")
    parser.add_option("-l", "--lower_bound", dest="lower",
		      help="lower bound for BSR comparisons, defaults to 0.4",
		      default="0.4", type="float")

    options, args = parser.parse_args()
    
    mandatories = ["matrix", "group1", "group2", "fasta"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.matrix,options.group1,options.group2,options.fasta,options.upper,options.lower)
