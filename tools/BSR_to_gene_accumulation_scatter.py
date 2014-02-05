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
        print "option not supported.  Only select acc, uni, or all"
        sys.exit()

def process_pangenome(matrix, upper, lower, iterations, type):
    my_matrix = open(matrix, "U")
    if type == "acc":
        acc_outfile = open("accumulation_replicates.txt", "w")
    elif type == "uni":
        uni_outfile = open("uniques_replicates.txt", "w")
    elif type == "core":
        core_outfile = open("core_replicates.txt", "w")
    else:
        acc_outfile = open("accumulation_replicates.txt", "w")
        uni_outfile = open("uniques_replicates.txt", "w")
        core_outfile = open("core_replicates.txt", "w")
    firstLine = my_matrix.readline()
    first_fields = firstLine.split()
    genomes = len(first_fields)
    indexes = []
    for x in first_fields:
        indexes.append(first_fields.index(x)+1)
    my_matrix.close()
    acc_dict = {}
    core_dict = {}
    uni_dict = {}
    for j in range(1,iterations+1):
        for i in range(1,genomes+1):
            positives_acc = []
            positives_core = []
            positives_unis = []
            outseqs=random.sample(set(indexes), int(i))
            with open(matrix, "U") as f:
                next(f)
                for line in f:
                    fields = line.split()
                    positive_lines_acc=[]
                    positive_lines_core=[]
                    positive_lines_unis=[]
                    for outseq in outseqs:
                        if type == "acc" or type == "all":
                            if float(fields[outseq])>=float(upper):
                                positive_lines_acc.append("1")
                        if type == "core" or type == "all":
                            if float(fields[outseq])>=float(upper):
                                positive_lines_core.append("1")
                        if type == "uni" or type == "all":
                            """this was changed from lower to upper"""
                            if float(fields[outseq])>=float(lower) and float(fields[outseq])>=float(upper):
                                positive_lines_unis.append("1")
                    if len(positive_lines_acc)>=1:
                        positives_acc.append("1")
                    if len(positive_lines_core)==len(outseqs):
                        positives_core.append("1")
                    if int(len(positive_lines_unis))==1:
                        positives_unis.append("1")
                    positive_lines_acc=[]
                    positive_lines_core=[]
                    positive_lines_unis=[]
            try:
                acc_dict[i].append(len(positives_acc))
            except KeyError:
                acc_dict[i] = [len(positives_acc)]
            try:
                core_dict[i].append(len(positives_core))
            except KeyError:
                core_dict[i] = [len(positives_core)]
            try:
                uni_dict[i].append(len(positives_unis))
            except KeyError:
                uni_dict[i] = [len(positives_unis)]
    try:
        sorted_acc_dict = collections.OrderedDict(sorted(acc_dict.items()))
        sorted_uni_dict = collections.OrderedDict(sorted(uni_dict.items()))
        sorted_core_dict = collections.OrderedDict(sorted(core_dict.items()))
    except:
        pass
    if type == "acc" or type == "all":
        print "accumulation means"
        for k,v in sorted_acc_dict.iteritems():
            print k, sum(v)/len(v)
            for z in v:
                print >> acc_outfile, str(k)+"\t"+str(z)+"\n",
    if type == "uni" or type == "all":
        print "unique means"
        for k,v in sorted_uni_dict.iteritems():
            print k, sum(v)/len(v)
            for z in v:
                print >> uni_outfile, str(k)+"\t"+str(z)+"\n",
    if type == "core" or type == "all":
        print "core means"
        for k,v in sorted_core_dict.iteritems():
            print k, sum(v)/len(v)
            for z in v:
                print >> core_outfile, str(k)+"\t"+str(z)+"\n",
    try:
        acc_outfile.close()
        uni_outfile.close()
        core_outfile.close()
    except:
        pass
    
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
