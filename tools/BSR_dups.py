#!/usr/bin/env python
"""Identifies duplicates in a BSR matrix"""

import sys
import optparse
import os
import errno
from ls_bsr.util import translate_genes
from ls_bsr.util import blast_against_self
from ls_bsr.util import parse_self_blast
from ls_bsr.util import blast_against_each_genome
from ls_bsr.util import get_seq_name
from igs.utils import logging
import subprocess
import glob

def filter_blast_report(ref_scores, length, max_plog, min_hlog):
    curr_dir=os.getcwd()
    orthologs = [ ]
    paralogs = [ ]
    for infile in glob.glob(os.path.join(curr_dir, "*_blast.out")):
        my_dict_o = {}
        my_dict_p = {}
        for line in open(infile, "U"):
            fields = line.split()
            if float(fields[2])>=float(min_hlog) and (float(fields[11])/float(ref_scores.get(fields[0])))>=float(length):
                my_dict_o.update({fields[0]:fields[11]})
        print my_dict_o
                #for k, v in my_dict:
                #for x in v:
                #print len(x)
            
def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print "%s cannot be found" % option
        sys.exit()

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

def main(directory, processors, genes, penalty, reward, filter, length, max_plog, min_hlog):
    start_dir = os.getcwd()
    ap=os.path.abspath("%s" % start_dir)
    dir_path=os.path.abspath("%s" % directory)
    try:
        os.makedirs('%s/joined' % dir_path)
    except OSError, e:
     	if e.errno != errno.EEXIST:
            raise
    for infile in glob.glob(os.path.join(dir_path, '*.fasta')):
        name=get_seq_name(infile)
        os.system("cp %s %s/joined/%s.new" % (infile,dir_path,name))
    gene_path=os.path.abspath("%s" % genes)
    os.system("cp %s %s/joined/" % (gene_path,dir_path))
    os.chdir("%s/joined" % dir_path)
    translate_genes(gene_path,)
    try:
        subprocess.check_call("formatdb -i %s -p F" % gene_path, shell=True)
    except:
        logging.logPrint("BLAST not found")
        sys.exit()
    blast_against_self(gene_path, "genes.pep", "tmp_blast.out", filter, "tblastn", penalty, reward)
    subprocess.check_call("sort -u -k 1,1 tmp_blast.out > self_blast.out", shell=True)
    ref_scores=parse_self_blast(open("self_blast.out", "U"))
    subprocess.check_call("rm tmp_blast.out self_blast.out", shell=True)
    logging.logPrint("starting BLAST")
    blast_against_each_genome(dir_path, processors, filter, "genes.pep", "tblastn", penalty, reward)
    filter_blast_report(ref_scores, length, max_plog, min_hlog)
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-d", "--directory", dest="directory",
                      help="/path/to/fasta_directory [REQUIRED]",
                      type="string", action="callback", callback=test_dir)
    parser.add_option("-p", "--parallel_workers", dest="processors",
                      help="How much work to do in parallel, defaults to 2, should number of CPUs your machine has",
                      default="2", type="int")
    parser.add_option("-g", "--genes", dest="genes", action="callback", callback=test_file,
                      help="predicted genes (nucleotide) to screen against genomes, will not use prodigal [REQUIRED]",
                      type="string")
    parser.add_option("-q", "--penalty", dest="penalty", action="store",
                      help="mismatch penalty, only to be used with blastn and -g option, default is -4",
                      default="-4", type="int")
    parser.add_option("-r", "--reward", dest="reward", action="store",
                      help="match reward, only to be used with blastn and -g option, default is 5",
                      default="5", type="int")
    parser.add_option("-f", "--filter", dest="filter", action="callback", callback=test_filter,
                      help="to use blast filtering or not, default is F or filter, change to T to turn off filtering",
                      default="F", type="string")
    parser.add_option("-l", "--length", dest="length", action="store",
                      help="Minimum length percentage of peptide to be called a homolog or paralog",
                      default="0.9", type="float")
    parser.add_option("-m", "--m_plog", dest="max_plog", action="store",
                      help="maximum value to be called a paralog, defaults to 0.80",
                      default="0.8", type="float")
    parser.add_option("-o", "--m_hlog", dest="min_hlog", action="store",
                      help="minimum value to be called a ortholog, defaults to 0.90",
                      default="0.9", type="float")
    options, args = parser.parse_args()
    
    mandatories = ["directory", "genes"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.directory, options.processors, options.genes, options.penalty,
         options.reward, options.filter, options.length, options.max_plog, options.min_hlog)
