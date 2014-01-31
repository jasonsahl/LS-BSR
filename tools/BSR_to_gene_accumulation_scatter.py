#!/usr/bin/env python

"""Sub-samples a LS-BSR matrix
and creates data that can generate
a scatterplot in Excel"""
from __future__ import division
from optparse import OptionParser
import sys
import random
import collections

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
        print "option not supported.  Only select acc, uni, or both"
        sys.exit()

def core_accumulation(matrix, upper, iterations):
    my_matrix = open(matrix, "U")
    outfile = open("core_replicates.txt", "w")
    firstLine = my_matrix.readline()
    first_fields = firstLine.split()
    genomes = len(first_fields)
    indexes = []
    for x in first_fields:
        indexes.append(first_fields.index(x)+1)
    my_matrix.close()
    total_dict = {}
    for j in range(0,iterations):
        for i in range(1,genomes+1):
            positives = [ ]
            outseqs=random.sample(set(indexes), int(i))
            with open(matrix, "U") as f:
                next(f)
                for line in f:
                    fields = line.split()
                    positive_lines=[]
                    for outseq in outseqs:
                        if float(fields[outseq])>=float(upper):
                             positive_lines.append("1")
                    if len(positive_lines)==len(outseqs):
                        positives.append("1")
            try:
                total_dict[i].append(len(positives))
            except KeyError:
                total_dict[i] = [len(positives)]
    sorted_dict = collections.OrderedDict(sorted(total_dict.items()))
    for k,v in sorted_dict.iteritems():
        for z in v:
            print >> outfile, str(k)+"\t"+str(z)+"\n",
    outfile.close()
       
def gene_accumulation(matrix, upper, iterations):
    my_matrix = open(matrix, "U")
    outfile = open("accumulation_replicates.txt", "w")
    firstLine = my_matrix.readline()
    first_fields = firstLine.split()
    genomes = len(first_fields)
    indexes = []
    for x in first_fields:
        indexes.append(first_fields.index(x)+1)
    my_matrix.close()
    total_dict = {}
    for j in range(0,iterations):
        for i in range(1,genomes+1):
            positives = [ ]
            outseqs=random.sample(set(indexes), int(i))
            with open(matrix, "U") as f:
                next(f)
                for line in f:
                    fields = line.split()
                    positive_lines=[]
                    for outseq in outseqs:
                        if float(fields[outseq])>=float(upper):
                             positive_lines.append("1")
                    if len(positive_lines)>=1:
                        positives.append("1")
            try:
                total_dict[i].append(len(positives))
            except KeyError:
                total_dict[i] = [len(positives)]
    sorted_dict = collections.OrderedDict(sorted(total_dict.items()))
    for k,v in sorted_dict.iteritems():
        for z in v:
            print >> outfile, str(k)+"\t"+str(z)+"\n",
    outfile.close()
    
def gene_uniques(matrix, upper, lower, iterations):
    my_matrix = open(matrix, "U")
    outfile = open("uniques_replicates.txt", "w")
    firstLine = my_matrix.readline()
    first_fields = firstLine.split()
    genomes = len(first_fields)
    indexes = []
    for x in first_fields:
        indexes.append(first_fields.index(x)+1)
    my_matrix.close()
    total_dict = {}
    for j in range(0,iterations):
        for i in range(1,genomes+1):
            positives = [ ]
            outseqs=random.sample(set(indexes), int(i))
            with open(matrix, "U") as f:
                next(f)
                for line in f:
                    fields = line.split()
                    positive_lines=[]
                    #for field in fields[1:]:
                    #    if float(field)>=float(lower):
                    #        positive_lines.append("1")
                    for outseq in outseqs:
                        if float(fields[outseq])>=float(lower):
                            positive_lines.append("1")
                    if len(positive_lines)==1:
                        positives.append("1")
            try:
                total_dict[i].append(len(positives))
            except KeyError:
                total_dict[i] = [len(positives)]
    sorted_dict = collections.OrderedDict(sorted(total_dict.items()))
    for k,v in sorted_dict.iteritems():
        for z in v:
            print >> outfile, str(k)+"\t"+str(z)+"\n",
    outfile.close()
    
def main(matrix, upper, lower, iterations, type):
    if "core" == type:
        core_accumulation(matrix, upper, iterations)
    elif "acc" == type:
        gene_accumulation(matrix, upper, iterations)
    elif "uni" == type:
        gene_uniques(matrix, upper, lower, iterations)
    else:
        core_accumulation(matrix, upper, iterations)
        gene_accumulation(matrix, upper, iterations)
        gene_uniques(matrix, upper, lower, iterations)
        
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
