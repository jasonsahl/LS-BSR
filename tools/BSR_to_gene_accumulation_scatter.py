#!/usr/bin/env python

"""Sub-samples a LS-BSR matrix
and creates data that can generate
a scatterplot in Excel"""
from optparse import OptionParser
import sys
from ls_bsr.util import process_pangenome

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def test_types(option, opt_str, value, parser):
    if "acc" in value:
        setattr(parser.values, option.dest, value)
    elif "uni" in value:
        setattr(parser.values, option.dest, value)
    elif "core" in value:
        setattr(parser.values, option.dest, value)
    elif "all" in value:
        setattr(parser.values, option.dest, value)
    else:
        print "option not supported.  Only select acc, uni, or all"
        sys.exit()

def main(matrix, upper, lower, iterations, type):
    process_pangenome(matrix, upper, lower, iterations, type)
            
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage) 
    parser.add_option("-b", "--bsr_matrix", dest="matrix",
                      help="path to BSR matrix [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-u", "--upper_bound", dest="upper",
                      help="upper bound to be called conserved, defaults to 0.8",
                      default="0.8", type="float")
    parser.add_option("-l", "--lower_bound", dest="lower",
                      help="lower bound to be called conserved, defaults to 0.4",
                      default="0.4", type="float")
    parser.add_option("-n", "--iterations", dest="iterations",
                      help="number of random samplings to perform, defaults to 10",
                      default="10", type="int", action="store")
    parser.add_option("-t", "--type", dest="type",
                      help="run accumulation (acc), uniques (uni), core(core), or all; defaults to all",
                      action="callback", callback=test_types, default="all", type="string")
    
    options, args = parser.parse_args()
    
    mandatories = ["matrix"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.matrix,options.upper,options.lower,options.iterations,options.type)
