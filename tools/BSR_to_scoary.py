"""BSR to Scoary"""

#!/usr/bin/env python

from __future__ import print_function
from optparse import OptionParser
import sys

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s file cannot be opened' % option)
        sys.exit()

def main(matrix,lower):
    my_matrix = open(matrix, "U")
    outfile = open("Scoary_matrix.txt","w")
    firstLine = my_matrix.readline()
    first_fields = firstLine.split()
    """Need to insert three blank fields, all comma separated"""
    first_fields.insert(0,"")
    first_fields.insert(1,"")
    first_fields.insert(2,"")
    outfile.write(",".join(first_fields)+"\n")
    for line in my_matrix:
        new_fields = []
        fields = line.split()
        new_fields.append(fields[0])
        """Adding in two blank fields"""
        new_fields.append("")
        new_fields.append("")
        for x in fields[1:]:
            if float(x)>=float(lower):
                new_fields.append("1")
            else:
                new_fields.append("")
        outfile.write(",".join(new_fields)+"\n")
    my_matrix.close()
    outfile.close()

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-b", "--bsr_matrix", dest="matrix",
                      help="path to BSR matrix [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-l", "--lower_bound", dest="lower",
                      help="lower bound to be called conserved, defaults to >=0.8",
                      default="0.8", type="float")
    options, args = parser.parse_args()

    mandatories = ["matrix"]
    for m in mandatories:
        if not options.__dict__[m]:
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)

    main(options.matrix,options.lower)
