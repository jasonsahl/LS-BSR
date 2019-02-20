#!/usr/bin/env python

"""extract only the unique IDs from a BSR matrix"""

from __future__ import print_function
from optparse import OptionParser
import sys
import os

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s file cannot be opened' % option)
        sys.exit()

def filter_uniques(matrix, threshold):
    outfile = open("uniques_BSR_matrix", "w")
    with open(matrix) as in_matrix:
        firstLine = in_matrix.readline()
    outdata = []
    with open(matrix) as in_matrix:
        for line in in_matrix:
            fields = line.split()
            totals = len(fields[1:])
            presents = []
            for x in fields[1:]:
                try:
                    if float(x)>=float(threshold):
                        presents.append(fields[0])
                except:
                    pass
            if int(len(presents))<int(2):
                outdata.append(fields[0])
                outfile.write(line)
    outfile.close()
    return outdata


def main(matrix, threshold):
    filter_uniques(matrix, threshold)

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-b", "--bsr_matrix", dest="matrix",
                      help="/path/to/bsr_matrix [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-t", "--threshold", dest="threshold",
                      help="lower threshold for ORF presence, defaults to 0.4",
                      action="store", default="0.4", type="float")

    options, args = parser.parse_args()

    mandatories = ["matrix"]
    for m in mandatories:
        if not options.__dict__[m]:
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)

    main(options.matrix, options.threshold)
