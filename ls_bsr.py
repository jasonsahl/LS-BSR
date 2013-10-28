#!/usr/bin/env python

"""Calculate the BSR value for all predicted CDSs
in a set of genomes in fasta format.

written by Jason Sahl
contacted at jasonsahl@gmail.com
"""

import sys
import os
import optparse
import subprocess
from subprocess import call
import errno
import types
from ls_bsr.util import *
from igs.utils import logging

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

def test_fplog(option, opt_str, value, parser):
    if "F" in value:
        setattr(parser.values, option.dest, value)
    elif "T" in value:
        setattr(parser.values, option.dest, value)
    else:
        print "select from T or F for f_plog setting"
        sys.exit()

def main(directory, id, filter, processors, genes, usearch, blast, penalty, reward, length,
         max_plog, min_hlog, f_plog, keep):
    start_dir = os.getcwd()
    ap=os.path.abspath("%s" % start_dir)
    dir_path=os.path.abspath("%s" % directory)
    ab = subprocess.call(['which', 'blastall'])
    if ab == 0:
        pass
    else:
        print "blastall isn't in your path, but needs to be!"
    try:
        os.makedirs('%s/joined' % dir_path)
    except OSError, e:
     	if e.errno != errno.EEXIST:
            raise
    for infile in glob.glob(os.path.join(dir_path, '*.fasta')):
        name=get_seq_name(infile)
        os.system("cp %s %s/joined/%s.new" % (infile,dir_path,name))
    if "null" in genes:
        rc = subprocess.call(['which', 'prodigal'])
        if rc == 0:
            pass
        else:
            print "prodigal is not in your path, but needs to be!"
        try:
            if os.path.exists(usearch):
                pass
        except:
            raise TypeError("-u usearch flag must be set for use with prodigal")
            sys.exc_clear()
        logging.logPrint("predicting genes with Prodigal")
        predict_genes(dir_path, processors)
        logging.logPrint("Prodigal done")
        os.system("cat *genes.seqs > all_gene_seqs.out")
        filter_scaffolds("all_gene_seqs.out")
        os.system("mv tmp.out all_gene_seqs.out")
        rename_fasta_header("all_gene_seqs.out", "all_sorted.txt")
        os.system("mkdir split_files")
        os.system("cp all_sorted.txt split_files/")
        os.system("rm all_sorted.txt")
        os.chdir("split_files/")
        os.system("split -l 200000 all_sorted.txt")
        sort_usearch(usearch)
        run_usearch(usearch, id)
        os.system("cat *.usearch.out > all_sorted.txt")
        os.system("mv all_sorted.txt %s/joined" % dir_path)
        os.chdir("%s/joined" % dir_path)
        uclust_cluster(usearch, id)
        translate_consensus("consensus.fasta")
        filter_seqs("tmp.pep")
        subprocess.check_call("formatdb -i consensus.fasta -p F", shell=True)
        blast_against_self("consensus.fasta", "consensus.pep", "tmp_blast.out", filter, blast, penalty, reward)
        subprocess.check_call("sort -u -k 1,1 tmp_blast.out > self_blast.out", shell=True)
        ref_scores=parse_self_blast(open("self_blast.out", "U"))
        subprocess.check_call("rm tmp_blast.out self_blast.out", shell=True)
        os.system("rm *new_genes.*")
        logging.logPrint("starting BLAST")
        blast_against_each_genome(dir_path, processors, filter, "consensus.pep", blast, penalty, reward)
        find_dups(ref_scores, length, max_plog, min_hlog)
    else:
        logging.logPrint("Using pre-compiled set of predicted genes")
        gene_path=os.path.abspath("%s" % genes)
        os.system("cp %s %s/joined/" % (gene_path,dir_path))
        os.chdir("%s/joined" % dir_path)
        if "tblastn" in blast:
            logging.logPrint("using tblastn")
            translate_genes(gene_path,)
            try:
                subprocess.check_call("formatdb -i %s -p F" % gene_path, shell=True)
            except:
                logging.logPrint("problem encountered with BLAST database")
                sys.exit()
            blast_against_self(gene_path, "genes.pep", "tmp_blast.out", filter, blast, penalty, reward)
            subprocess.check_call("sort -u -k 1,1 tmp_blast.out > self_blast.out", shell=True)
            ref_scores=parse_self_blast(open("self_blast.out", "U"))
            subprocess.check_call("rm tmp_blast.out self_blast.out", shell=True)
            logging.logPrint("starting BLAST")
            blast_against_each_genome(dir_path, processors, filter, "genes.pep", blast, penalty, reward)
        else:
            logging.logPrint("using blastn")
            try:
                subprocess.check_call("formatdb -i %s -p F" % gene_path, shell=True)
            except:
                logging.logPrint("BLAST not found")
                sys.exit()
            blast_against_self(gene_path, gene_path, "tmp_blast.out", filter, blast, penalty, reward)
            subprocess.check_call("sort -u -k 1,1 tmp_blast.out > self_blast.out", shell=True)
            ref_scores=parse_self_blast(open("self_blast.out", "U"))
            subprocess.check_call("rm tmp_blast.out self_blast.out", shell=True)
            logging.logPrint("starting BLAST")
            blast_against_each_genome(dir_path, processors, filter, gene_path, blast, penalty, reward)
    logging.logPrint("BLAST done")
    parse_blast_report()
    get_unique_lines()
    make_table(processors, "F")
    subprocess.check_call("paste ref.list *.matrix > bsr_matrix", shell=True)
    divide_values("bsr_matrix", ref_scores)
    subprocess.check_call("paste ref.list BSR_matrix_values.txt > %s/bsr_matrix_values.txt" % start_dir, shell=True)
    if "T" in f_plog:
        filter_paralogs("%s/bsr_matrix_values.txt" % start_dir, "paralog_ids.txt")
        os.system("cp bsr_matrix_values_filtered.txt %s" % start_dir)
    else:
        pass
    try:
        subprocess.check_call("cp names.txt consensus.pep consensus.fasta duplicate_ids.txt paralog_ids.txt %s" % start_dir, shell=True, stderr=open(os.devnull, 'w'))
    except:
        sys.exc_clear()
    logging.logPrint("all Done")
    os.chdir("%s" % dir_path)
    if "T" == keep:
        pass
    else:
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
                      help="How much work to do in parallel, defaults to 2",
                      default="2", type="int")
    parser.add_option("-g", "--genes", dest="genes", action="callback", callback=test_file,
                      help="predicted genes (nucleotide) to screen against genomes, will not use prodigal",
                      type="string",default="null")
    parser.add_option("-u", "--usearch", dest="usearch", action="store",
                      help="path to usearch v6, required for use with Prodigal",
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
    parser.add_option("-l", "--length", dest="length", action="store",
                      help="minimum BSR value to be called a duplicate, defaults to 0.7",
                      default="0.7", type="float")
    parser.add_option("-m", "--max_plog", dest="max_plog", action="store",
                      help="maximum value to be called a paralog, defaults to 0.85",
                      default="0.85", type="float")
    parser.add_option("-n", "--min_hlog", dest="min_hlog", action="store",
                      help="minimum BLAST ID to be called a homolog, defaults to 75",
                      default="75", type="int")
    parser.add_option("-t", "--f_plog", dest="f_plog", action="callback", callback=test_fplog,
                      help="filter ORFs with a paralog from BSR matrix? Default is F, values can be T or F",
                      default="F", type="string")
    parser.add_option("-k", "--keep", dest="keep", action="callback", callback=test_filter,
                      help="keep or remove temp files, choose from T or F, defaults to F",
                      default="F", type="string")
    options, args = parser.parse_args()
    
    mandatories = ["directory"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.directory, options.id, options.filter, options.processors, options.genes, options.usearch, options.blast,
         options.penalty, options.reward, options.length, options.max_plog, options.min_hlog, options.f_plog, options.keep)

