#!/usr/bin/env python

"""removes conserve ORFs, so we are only looking at
the variable region of the pan-genome"""

from optparse import OptionParser
from ls_bsr.util import filter_variome

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def main(matrix, threshold, step):
    filter_variome(matrix, threshold, step)

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-b", "--bsr_matrix", dest="matrix",
                      help="/path/to/bsr_matrix [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-t", "--threshold", dest="threshold",
                      help="lower threshold for ORF presence, defaults to 0.8",
                      action="store", default="0.8", type="float")
    parser.add_option("-s", "--step", dest="step",
                      help="how many genomes fewer than the total to consider for loss, defaults to 1",
                      action="store", default="1", type="int")
   
    options, args = parser.parse_args()
    
    mandatories = ["matrix"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.matrix, options.threshold, options.step)
