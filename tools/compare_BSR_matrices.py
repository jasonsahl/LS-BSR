#!/usr/bin/env python

"""compares BSR values between two groups in a BSR matrix
Numpy and BioPython need to be installed.  Python version must be at
least 2.7 to use collections"""

from optparse import OptionParser
from collections import deque
import numpy as np
import subprocess
from Bio import SeqIO

def prune_matrix(matrix, group1, group2):
    """prune out genomes of interest from a BSR matrix.
    Not done efficiently, but appears to work"""
    in_matrix = open(matrix, "rU")
    group1_ids = [ ]
    group2_ids = [ ]
    group1_out = open("group1_pruned.txt", "w")
    group2_out = open("group2_pruned.txt", "w")
    for line in open(group1, "rU"):
        line.strip()
        group1_ids.append(line)
    for line in open(group2, "rU"):
        line.strip()
        group2_ids.append(line)
    firstLine = in_matrix.readline()
    fields = firstLine.split()
    fields.insert(0, "cluster")
    group1_ids = map(lambda s: s.strip(), group1_ids)
    group2_ids = map(lambda s: s.strip(), group2_ids)
    group1_idx = [ ]
    for x in fields:
        if x not in group1_ids: group1_idx.append(fields.index(x))  
    deque((list.pop(fields, i) for i in sorted(group1_idx, reverse=True)), maxlen=0)
    print >> group1_out, "\t","\t","\t".join(fields)
    for line in in_matrix:
        fields = line.split()
        name = fields[0]
        deque((list.pop(fields, i) for i in sorted(group1_idx, reverse=True)), maxlen=0)
	print >> group1_out,"".join(name),"\t","\t".join(fields)
    in_matrix = open(matrix, "rU")
    firstLine = in_matrix.readline()
    fields = firstLine.split()
    fields.insert(0, "cluster")
    group2_idx = [ ]
    for x in fields:
        if x not in group2_ids: group2_idx.append(fields.index(x))
    deque((list.pop(fields, i) for i in sorted(group2_idx, reverse=True)), maxlen=0)
    print >> group2_out, "\t", "\t".join(fields)
    for line in in_matrix:
        fields = line.split()
        name = fields[0]
        deque((list.pop(fields, i) for i in sorted(group2_idx, reverse=True)), maxlen=0)
        print >> group2_out,"".join(name),"\t","\t".join(fields)

def compare_values(pruned_1,pruned_2,upper,lower):
    group1 = open(pruned_1, "rU")
    group2 = open(pruned_2, "rU")
    group1_out = open("group1_out.txt", "w")
    group2_out = open("group2_out.txt", "w")
    next(group1)
    for line in group1:
	fields = line.split()
	presents = [ ]
	homolog = [ ]
	ints=map(float, fields[1:])
	mean = float(np.mean(ints))
	for x in fields[1:]:
		if float(x)>=float(upper): presents.append(x)
	for x in fields[1:]:
		if float(x)>=float(lower): homolog.append(x)
	print >> group1_out,fields[0],"\t",mean,"\t",len(presents),"\t",len(fields[1:]),"\t",len(homolog)
    next(group2)
    for line in group2:
	fields = line.split()
	presents = [ ]
	homolog = [ ]
	ints=map(float, fields[1:])
	mean = float(np.mean(ints))
	for x in fields[1:]:
		if float(x)>=float(upper): presents.append(x)
	for x in fields[1:]:
		if float(x)>=float(lower): homolog.append(x)
	print >> group2_out,mean,"\t",len(presents),"\t",len(fields[1:]),"\t",len(homolog)

def find_uniques(combined,fasta):
    infile = open(combined, "rU")
    group1_unique_ids = [ ]
    seqrecords=[ ]
    for line in infile:
	fields=line.split()
	if int(fields[2])/int(fields[3])==1 and int(fields[8])==0:
		group1_unique_ids.append(fields[0])
    for record in SeqIO.parse(fasta, "fasta"):
	    if record.id in group1_unique_ids:
		    seqrecords.append(record)
    output_handle = open("group1_unique_seqs.fasta", "w")
    SeqIO.write(seqrecords, output_handle, "fasta")
    output_handle.close()
    group2_unique_ids = [ ]
    seqrecords2 = [ ]
    infile = open(combined, "rU")
    for line in infile:
	fields=line.split()
	if int(fields[6])/int(fields[7])==1 and int(fields[4])==0:
	       group2_unique_ids.append(fields[0])
    for record in SeqIO.parse(fasta, "fasta"):
	    if record.id in group2_unique_ids:
		    seqrecords2.append(record)
    output_handle2 = open("group2_unique_seqs.fasta", "w")
    SeqIO.write(seqrecords2, output_handle2, "fasta")
    output_handle2.close()

        	
def main(matrix,group1,group2,fasta,upper,lower):
    prune_matrix(matrix,group1,group2)
    compare_values("group1_pruned.txt","group2_pruned.txt",upper,lower)
    subprocess.check_call("paste group1_out.txt group2_out.txt > groups_combined.txt", shell=True)
    find_uniques("groups_combined.txt",fasta)
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage) 
    parser.add_option("-b", "--bsr_matrix", dest="matrix",
                      help="/path/to/bsr_matrix [REQUIRED]",
                      action="store", type="string")
    parser.add_option("-f", "--fasta", dest="fasta",
                      help="/path/to/ORF_fasta_file [REQUIRED]",
                      action="store", type="string")
    parser.add_option("-1", "--group_1_ids", dest="group1",
                      help="new line separated file with group1 ids [REQUIRED]",
                      action="store", type="string")
    parser.add_option("-2", "--group_2_ids", dest="group2",
                      help="new line separated file with group2 ids [REQUIRED]",
                      action="store", type="string")
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
