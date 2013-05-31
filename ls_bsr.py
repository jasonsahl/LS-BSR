#!/usr/bin/env python

"""Calculate the BSR value for all predicted ORFs
in a set of genomes in fasta format.  V3 - replaced
transeq with BioPython.  V4 - changed to true BSR
test.  V5 - fixed bug in how BSR was calculated.
V6 - changed gene caller from Glimmer to Prodigal"""

import sys
import os
import optparse
import subprocess
from subprocess import call
import errno
import types
from ls_bsr.util import *


def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print 'genes file cannot be opened'
        sys.exit()

def test_filter(option, opt_str, value, parser):
    if "F" in value:
        setattr(parser.values, option.dest, value)
    elif "T" in value:
        setattr(parser.values, option.dest, value)
    else:
        print "option not supported.  Only select from T and F"
        sys.exit()

def test_blast(option, opt_str, value, parser):
    if "tblastn" in value:
        setattr(parser.values, option.dest, value)
    elif "blastn" in value:
        setattr(parser.values, option.dest, value)
    else:
        print "Blast option not supported.  Only select from blastn and blastn"
        sys.exit()

def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print "directory of fastas cannot be found"
        sys.exit()

def test_id(option, opt_str, value, parser):
    if type(value) == types.IntType:
        sys.exit()
    elif type(value) == types.FloatType:
        setattr(parser.values, option.dest, value)
    else:
        print "id value needs to be a float"
        sys.exit()

def test_usearch(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print "usearch can't be found"
        sys.exit()

def main(directory, id, filter, processors, genes, usearch, blast, penalty, reward):
    start_dir = os.getcwd()
    ap=os.path.abspath("%s" % directory)
    try:
     	os.makedirs('%s/joined' % directory)
    except OSError, e:
     	if e.errno != errno.EEXIST:
            raise
    for infile in glob.glob(os.path.join(directory, '*.fasta')):
        name=get_seq_name(infile)
        os.system("cp %s %s/joined/%s.new" % (infile,directory,name))
    if "null" in genes:
        logging.logPrint("predicting genes with Prodigal")
        predict_genes(directory, processors)
        logging.logPrint("Prodigal done")
        os.system("cat *genes.seqs > all_gene_seqs.out")
        uclust_sort(usearch)
        rename_fasta_header("tmp_sorted.txt", "all_sorted.txt")
        uclust_cluster(usearch, id)
        translate_consensus("consensus.fasta")
        filter_seqs("tmp.pep")
        subprocess.check_call("formatdb -i consensus.fasta -p F", shell=True)
        blast_against_self("consensus.fasta", "consensus.pep", "tmp_blast.out", filter, blast, penalty, reward)
        subprocess.check_call("sort -u -k 1,1 tmp_blast.out > self_blast.out", shell=True)
        ref_scores=parse_self_blast("self_blast.out")
        subprocess.check_call("rm tmp_blast.out self_blast.out", shell=True)
        os.system("rm *new_genes.*")
        logging.logPrint("starting BLAST")
        blast_against_each_genome(directory, processors, filter, "consensus.pep", blast, penalty, reward)
    else:
        logging.logPrint("Using pre-compiled set of predicted genes")
        os.system("cp %s %s/joined/" % (genes,directory))
        os.chdir("%s/joined" % directory)
        if "tblastn" in blast:
            logging.logPrint("using tblastn")
            translate_genes(genes, directory)
            subprocess.check_call("formatdb -i %s -p F" % genes, shell=True)
            blast_against_self(genes, "genes.pep", "tmp_blast.out", filter, blast, penalty, reward)
            subprocess.check_call("sort -u -k 1,1 tmp_blast.out > self_blast.out", shell=True)
            ref_scores=parse_self_blast("self_blast.out")
            subprocess.check_call("rm tmp_blast.out self_blast.out", shell=True)
            logging.logPrint("starting BLAST")
            blast_against_each_genome(directory, processors, filter, "genes.pep", blast, penalty, reward)
        else:
            logging.logPrint("using blastn")
            subprocess.check_call("formatdb -i %s -p F" % genes, shell=True)
            blast_against_self(genes, genes, "tmp_blast.out", filter, blast, penalty, reward)
            subprocess.check_call("sort -u -k 1,1 tmp_blast.out > self_blast.out", shell=True)
            ref_scores=parse_self_blast("self_blast.out")
            subprocess.check_call("rm tmp_blast.out self_blast.out", shell=True)
            logging.logPrint("starting BLAST")
            blast_against_each_genome(directory, processors, filter, genes, blast, penalty, reward)
    logging.logPrint("BLAST done")
    parse_blast_report(directory)
    get_unique_lines(directory)
    make_table(directory, processors)
    subprocess.check_call("paste ref.list *.matrix > bsr_matrix", shell=True)
    divide_values("bsr_matrix", ref_scores)
    subprocess.check_call("paste ref.list BSR_matrix_values.txt > %s/bsr_matrix_values.txt" % start_dir, shell=True)
    try:
        subprocess.check_call("cp names.txt consensus.pep consensus.fasta %s" % start_dir, shell=True, stderr=open(os.devnull, 'w'))
    except:
        sys.exc_clear()
    logging.logPrint("all Done")
    os.chdir("%s" % ap)
    os.system("rm -rf joined")
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-d", "--directory", dest="directory",
                      help="/path/to/fasta_directory [REQUIRED]",
                      type="string", action="callback", callback=test_dir)
    parser.add_option("-i", "--identity", dest="id", action="callback", callback=test_id,
                      help="clustering id for USEARCH (0.0-1.0), defaults to 0.9",
                      type="float", default="0.9")
    parser.add_option("-f", "--filter", dest="filter", action="callback", callback=test_filter,
                      help="to use blast filtering or not, default is F or filter, change to T to turn off filtering",
                      default="F", type="string")
    parser.add_option("-p", "--parallel_workers", dest="processors",
                      help="How much work to do in parallel, defaults to 2, should number of CPUs your machine has",
                      default="2", type="int")
    parser.add_option("-g", "--genes", dest="genes", action="callback", callback=test_file,
                      help="predicted genes (nucleotide) to screen against genomes, will not use prodigal",
                      type="string",default="null")
    parser.add_option("-u", "--usearch", dest="usearch", action="callback", callback=test_usearch,
                      help="path to usearch v6 [REQUIRED]",
                      type="string")
    parser.add_option("-b", "--blast", dest="blast", action="callback", callback=test_blast,
                      help="use tblast or blastn, only used in conjunction with -g option, default is blastn",
                      default="tblastn", type="string")
    parser.add_option("-q", "--penalty", dest="penalty", action="store",
                      help="mismatch penalty, only to be used with blastn and -g option, default is -4",
                      default="-4", type="int")
    parser.add_option("-r", "--reward", dest="reward", action="store",
                      help="match reward, only to be used with blastn and -g option, default is 5",
                      default="5", type="int")
                      
    options, args = parser.parse_args()
    
    mandatories = ["directory", "usearch"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.directory, options.id, options.filter, options.processors, options.genes, options.usearch, options.blast,
         options.penalty, options.reward)

