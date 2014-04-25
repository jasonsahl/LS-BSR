#!/usr/bin/env python

"""splits a reference sequence into chunks of N size,
with a provided step size"""

import itertools
from Bio import SeqIO
import os
import sys
import optparse

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s cannot be opened' % option
        sys.exit()

class Increments:
    def __init__(self, start, increment):
        self.state = start
        self.p_start = start
        self.p_increment = increment

    def next(self):
        self.state += self.p_increment
        return self.state

    def reset(self):
        self.state = self.p_start
        
record_count_1 = Increments(1, 1)
        
def get_seq_name(in_fasta):
    """used for renaming the sequences"""
    return os.path.basename(in_fasta)

def split_sequence_by_window(input_file, step_size, frag_length):
    """cuts up fasta sequences into given chunks"""
    infile = open(input_file, "rU")
    first_record = list(itertools.islice(SeqIO.parse(infile,"fasta"), 1))[0]
    return sliding_window(first_record.seq, frag_length, step_size)

def write_sequences(reads,name):
    """write shredded fasta sequences to disk"""
    handle = open("%s_seqs_shredded.txt" % name, "w")
    for read in reads:
        print >> handle, ">%d\n%s" % (record_count_1.next(), read)
    handle.close()

def sliding_window(sequence, frag_length, step_size=5):
    """cuts up sequence into a given length"""
    numOfChunks = (len(sequence) - frag_length) + 1
    for i in range(0, numOfChunks, step_size):
        yield sequence[i:i + frag_length]

def main(reference,frag_length,step_size):
    name = get_seq_name(reference)
    redux_name = name.replace(".fasta","")
    reads = split_sequence_by_window(reference, step_size, frag_length)
    write_sequences(reads, redux_name)
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-r", "--reference", dest="reference",
                      help="/path/to/reference_genome [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-f", "--frag_length", dest="frag_length",
                      help="length to slice fragment, defaults to 1000",
                      action="store", type="int", default="1000")
    parser.add_option("-s", "--step", dest="step_size",
                      help="step size for shredding sequences, defaults to 100",
                      default="100", type="int")
    options, args = parser.parse_args()

    mandatories = ["reference"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.reference,options.frag_length,options.step_size)

