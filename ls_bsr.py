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
import glob
import tempfile

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def test_filter(option, opt_str, value, parser):
    if "F" in value:
        setattr(parser.values, option.dest, value)
    elif "T" in value:
        setattr(parser.values, option.dest, value)
    else:
        print "option not supported.  Only select from T and F"
        sys.exit()

def test_cluster(option, opt_str, value, parser):
    if "usearch" in value:
        setattr(parser.values, option.dest, value)
    elif "vsearch" in value:
        setattr(parser.values, option.dest, value)
    elif "cd-hit" in value:
        setattr(parser.values, option.dest, value)
    else:
        print "option not supported. Choose from vsearch, usearch, cd-hit"
        sys.exit()

def test_blast(option, opt_str, value, parser):
    if "tblastn" in value:
        setattr(parser.values, option.dest, value)
    elif "blastn" in value:
        setattr(parser.values, option.dest, value)
    elif "blat" in value:
        setattr(parser.values, option.dest, value)
    elif "blastall" in value:
        setattr(parser.values, option.dest, value)
    else:
        print "Blast option not supported.  Only select from tblastn, blat, or blastn"
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

def test_fplog(option, opt_str, value, parser):
    if "F" in value:
        setattr(parser.values, option.dest, value)
    elif "T" in value:
        setattr(parser.values, option.dest, value)
    else:
        print "select from T or F for f_plog setting"
        sys.exit()

def main(directory,id,filter,processors,genes,cluster_method,blast,length,
         max_plog,min_hlog,f_plog,keep,filter_peps,filter_scaffolds,prefix,temp_dir,debug):
    start_dir = os.getcwd()
    ap=os.path.abspath("%s" % start_dir)
    dir_path=os.path.abspath("%s" % directory)
    logging.logPrint("Testing paths of dependencies")
    if blast=="blastn" or blast=="tblastn":
        ab = subprocess.call(['which', 'blastn'])
        if ab == 0:
            print "citation: Altschul SF, Madden TL, Schaffer AA, Zhang J, Zhang Z, Miller W, and Lipman DJ. 1997. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Nucleic Acids Res 25:3389-3402"
        else:
            print "blastn isn't in your path, but needs to be!"
            sys.exit()
    if "NULL" in temp_dir:
        fastadir = tempfile.mkdtemp()
    else:
        fastadir = os.path.abspath("%s" % temp_dir)
        if os.path.exists('%s' % temp_dir):
            print "old run directory exists in your genomes directory (%s).  Delete and run again" % temp_dir
            sys.exit()
        else:
            os.makedirs('%s' % temp_dir)
    for infile in glob.glob(os.path.join(dir_path, '*.fasta')):
        name=get_seq_name(infile)
        os.link("%s" % infile, "%s/%s.new" % (fastadir,name))
    if "null" in genes:
        rc = subprocess.call(['which', 'prodigal'])
        if rc == 0:
            pass
        else:
            print "prodigal is not in your path, but needs to be!"
            sys.exit()
        print "citation: Hyatt D, Chen GL, Locascio PF, Land ML, Larimer FW, and Hauser LJ. 2010. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11:119"
        if "usearch" in cluster_method:
            print "citation: Edgar RC. 2010. Search and clustering orders of magnitude faster than BLAST. Bioinformatics 26:2460-2461"
        elif "cd-hit" in cluster_method:
            print "citation: Li, W., Godzik, A. 2006. Cd-hit: a fast program for clustering and comparing large sets of protein or nuceltodie sequences. Bioinformatics 22(13):1658-1659"
        else:
            pass
        if blast=="blat":
            ac = subprocess.call(['which', 'blat'])
            if ac == 0:
                print "citation: W.James Kent. 2002. BLAT - The BLAST-Like Alignment Tool.  Genome Research 12:656-664"
            else:
                print "You have requested blat, but it is not in your PATH"
                sys.exit()
        logging.logPrint("predicting genes with Prodigal")
        #predict_genes(dir_path, processors)
        predict_genes(fastadir, processors)
        logging.logPrint("Prodigal done")
        """This function produces locus tags"""
        genbank_hits = process_genbank_files(dir_path)
        if genbank_hits == None or len(genbank_hits) == 0:
            os.system("cat *genes.seqs > all_gene_seqs.out")
            if filter_scaffolds == "T":
                filter_scaffolds("all_gene_seqs.out")
                os.system("mv tmp.out all_gene_seqs.out")
            else:
                pass
        else:
            logging.logPrint("Converting genbank files")
            """First combine all of the prodigal files into one file"""
            os.system("cat *genes.seqs > all_gene_seqs.out")
            if filter_scaffolds == "T":
                filter_scaffolds("all_gene_seqs.out")
                os.system("mv tmp.out all_gene_seqs.out")
            else:
                pass
            """This combines the locus tags with the Prodigal prediction"""
            os.system("cat *locus_tags.fasta all_gene_seqs.out > tmp.out")
            os.system("mv tmp.out all_gene_seqs.out")
            """I also need to convert the GenBank file to a FASTA file"""
            for hit in genbank_hits:
                reduced_hit = hit.replace(".gbk","")
                SeqIO.convert("%s/%s" % (dir_path, hit), "genbank", "%s.fasta.new" % reduced_hit, "fasta")
        if "NULL" in cluster_method:
            print "Clustering chosen, but no method selected...exiting"
            sys.exit()
        elif "usearch" in cluster_method:
            ac = subprocess.call(['which', 'usearch'])
            if ac == 0:
                os.system("mkdir split_files")
                os.system("mv all_gene_seqs.out split_files/all_sorted.txt")
                os.chdir("split_files/")
                logging.logPrint("Splitting FASTA file for use with USEARCH")
                split_files("all_sorted.txt")
                logging.logPrint("clustering with USEARCH at an ID of %s" % id)
                run_usearch(id)
                os.system("cat *.usearch.out > all_sorted.txt")
                #os.system("mv all_sorted.txt %s/joined" % dir_path)
                os.system("mv all_sorted.txt %s" % fastadir)
                os.chdir("%s" % fastadir)
                #os.chdir("%s/joined" % dir_path)
                uclust_cluster(id)
                logging.logPrint("USEARCH clustering finished")
            else:
                print "usearch must be in your path as usearch...exiting"
                sys.exit()
        elif "vsearch" in cluster_method:
            ac = subprocess.call(['which', 'vsearch'])
            if ac == 0:
                logging.logPrint("clustering with VSEARCH at an ID of %s, using %s processors" % (id,processors))
                run_vsearch(id, processors)
                os.system("mv vsearch.out consensus.fasta")
                logging.logPrint("VSEARCH clustering finished")
            else:
                print "vsearch must be in your path as vsearch...exiting"
                sys.exit()
        elif "cd-hit" in cluster_method:
            ac = subprocess.call(['which', 'cd-hit-est'])
            if ac == 0:
                logging.logPrint("clustering with cd-hit at an ID of %s, using %s processors" % (id,processors))
                subprocess.check_call("cd-hit-est -i all_gene_seqs.out -o consensus.fasta -M 0 -T %s -c %s > /dev/null 2>&1" % (processors, id), shell=True)
            else:
                print "cd-hit must be in your path as cd-hit-est...exiting"
                sys.exit()
        """need to check for dups here"""
        dup_ids = test_duplicate_header_ids("consensus.fasta")
        if dup_ids == "True":
            pass
        elif dup_ids == "False":
            print "duplicate headers identified, renaming.."
            rename_fasta_header("consensus.fasta", "tmp.txt")
            os.system("mv tmp.txt consensus.fasta")
        else:
            pass
        if "tblastn" == blast:
            subprocess.check_call("makeblastdb -in consensus.fasta -dbtype nucl > /dev/null 2>&1", shell=True)
            translate_consensus("consensus.fasta")
            if filter_peps == "T":
                filter_seqs("tmp.pep")
                os.system("rm tmp.pep")
            else:
                os.system("mv tmp.pep consensus.pep")
            clusters = get_cluster_ids("consensus.pep")
            blast_against_self_tblastn("tblastn", "consensus.fasta", "consensus.pep", "tmp_blast.out", processors, filter)
        elif "blastn" == blast:
            subprocess.check_call("makeblastdb -in consensus.fasta -dbtype nucl > /dev/null 2>&1", shell=True)
            blast_against_self_blastn("blastn", "consensus.fasta", "consensus.fasta", "tmp_blast.out", filter, processors)
            clusters = get_cluster_ids("consensus.fasta")
        elif "blat" == blast:
            blat_against_self("consensus.fasta", "consensus.fasta", "tmp_blast.out", processors)
            clusters = get_cluster_ids("consensus.fasta")
        else:
            pass
        subprocess.check_call("sort -u -k 1,1 tmp_blast.out > self_blast.out", shell=True)
        ref_scores=parse_self_blast(open("self_blast.out", "U"))
        subprocess.check_call("rm tmp_blast.out self_blast.out", shell=True)
        os.system("rm *new_genes.*")
        if blast == "tblastn" or blast == "blastn":
            logging.logPrint("starting BLAST")
        else:
            logging.logPrint("starting BLAT")
        if "tblastn" == blast:
            blast_against_each_genome_tblastn(dir_path, processors, "consensus.pep", filter)
        elif "blastn" == blast:
            blast_against_each_genome_blastn(dir_path, processors, filter, "consensus.fasta")
        elif "blat" == blast:
            blat_against_each_genome(dir_path, "consensus.fasta",processors)
        else:
            pass
    else:
        logging.logPrint("Using pre-compiled set of predicted genes")
        files = glob.glob(os.path.join(dir_path, "*.fasta"))
        if len(files)==0:
            print "no usable reference genomes found!"
            sys.exit()
        else:
            pass
        gene_path=os.path.abspath("%s" % genes)
        dup_ids = test_duplicate_header_ids(gene_path)
        if dup_ids == "True":
            pass
        elif dup_ids == "False":
            print "duplicate headers identified, exiting.."
            sys.exit()
        clusters = get_cluster_ids(gene_path)
        #os.system("cp %s %s/joined/" % (gene_path,dir_path))
        os.system("cp %s %s" % (gene_path,fastadir))
        #os.chdir("%s/joined" % dir_path)
        os.chdir("%s" % fastadir)
        if gene_path.endswith(".pep"):
            logging.logPrint("using tblastn on peptides")
            try:
                subprocess.check_call("makeblastdb -in %s -dbtype prot > /dev/null 2>&1" % gene_path, shell=True)
            except:
                logging.logPrint("problem encountered with BLAST database")
                sys.exit()
            blast_against_self_tblastn("blastp", gene_path, gene_path, "tmp_blast.out", processors, filter)
            subprocess.check_call("sort -u -k 1,1 tmp_blast.out > self_blast.out", shell=True)
            ref_scores=parse_self_blast(open("self_blast.out", "U"))
            subprocess.check_call("rm tmp_blast.out self_blast.out", shell=True)
            logging.logPrint("starting BLAST")
            blast_against_each_genome_tblastn(dir_path, processors, gene_path, filter)
        elif gene_path.endswith(".fasta"):
            if "tblastn" == blast:
                logging.logPrint("using tblastn")
                translate_genes(gene_path)
                try:
                    subprocess.check_call("makeblastdb -in %s -dbtype nucl > /dev/null 2>&1" % gene_path, shell=True)
                except:
                    logging.logPrint("problem encountered with BLAST database")
                    sys.exit()
                blast_against_self_tblastn("tblastn", gene_path, "genes.pep", "tmp_blast.out", processors, filter)
                subprocess.check_call("sort -u -k 1,1 tmp_blast.out > self_blast.out", shell=True)
                ref_scores=parse_self_blast(open("self_blast.out", "U"))
                subprocess.check_call("rm tmp_blast.out self_blast.out", shell=True)
                logging.logPrint("starting BLAST")
                blast_against_each_genome_tblastn(dir_path, processors, "genes.pep", filter)
                os.system("cp genes.pep %s" % start_dir)
            elif "blastn" == blast:
                logging.logPrint("using blastn")
                try:
                    subprocess.check_call("makeblastdb -in %s -dbtype nucl > /dev/null 2>&1" % gene_path, shell=True)
                except:
                    logging.logPrint("Database not formatted correctly...exiting")
                    sys.exit()
                try:
                    blast_against_self_blastn("blastn", gene_path, gene_path, "tmp_blast.out", filter, processors)
                except:
                    print "problem with blastn, exiting"
                    sys.exit()
                subprocess.check_call("sort -u -k 1,1 tmp_blast.out > self_blast.out", shell=True)
                ref_scores=parse_self_blast(open("self_blast.out", "U"))
                subprocess.check_call("rm tmp_blast.out self_blast.out", shell=True)
                logging.logPrint("starting BLAST")
                try:
                    blast_against_each_genome_blastn(dir_path, processors, filter, gene_path)
                except:
                    print "problem with blastn, exiting"
                    sys.exit()
            elif "blat" == blast:
                logging.logPrint("using blat")
                blat_against_self(gene_path, gene_path, "tmp_blast.out", processors)
                subprocess.check_call("sort -u -k 1,1 tmp_blast.out > self_blast.out", shell=True)
                ref_scores=parse_self_blast(open("self_blast.out", "U"))
                subprocess.check_call("rm tmp_blast.out self_blast.out", shell=True)
                logging.logPrint("starting BLAT")
                blat_against_each_genome(dir_path,gene_path,processors)
            else:
                pass
        else:
            print "input file format not supported"
            sys.exit()
    find_dups_dev(ref_scores, length, max_plog, min_hlog, clusters, processors)
    if blast=="blat":
        logging.logPrint("BLAT done")
    else:
        logging.logPrint("BLAST done")
    parse_blast_report("false")
    get_unique_lines()
    curr_dir=os.getcwd()
    table_files = glob.glob(os.path.join(curr_dir, "*.filtered.unique"))
    files_and_temp_names = [(str(idx), os.path.join(curr_dir, f))
                            for idx, f in enumerate(table_files)]
    names=[]
    table_list = []
    nr_sorted=sorted(clusters)
    centroid_list = []
    centroid_list.append(" ")
    for x in nr_sorted:
        centroid_list.append(x)
    table_list.append(centroid_list)
    logging.logPrint("starting matrix building")
    new_names,new_table = new_loop(files_and_temp_names, processors, clusters, debug)
    new_table_list = table_list+new_table
    logging.logPrint("matrix built")
    open("ref.list", "a").write("\n")
    for x in nr_sorted:
        open("ref.list", "a").write("%s\n" % x)
    names_out = open("names.txt", "w")
    names_redux = [val for subl in new_names for val in subl]
    for x in names_redux: print >> names_out, "".join(x)
    names_out.close()
    create_bsr_matrix_dev(new_table_list)
    divide_values("bsr_matrix", ref_scores)
    subprocess.check_call("paste ref.list BSR_matrix_values.txt > %s/bsr_matrix_values.txt" % start_dir, shell=True)
    if "T" in f_plog:
        filter_paralogs("%s/bsr_matrix_values.txt" % start_dir, "paralog_ids.txt")
        os.system("cp bsr_matrix_values_filtered.txt %s" % start_dir)
    else:
        pass
    try:
        subprocess.check_call("cp dup_matrix.txt names.txt consensus.pep consensus.fasta duplicate_ids.txt paralog_ids.txt %s" % ap, shell=True, stderr=open(os.devnull, 'w'))
    except:
        sys.exc_clear()
    """new code to rename files according to a prefix"""
    import datetime
    timestamp = datetime.datetime.now()
    rename = str(timestamp.year), str(timestamp.month), str(timestamp.day), str(timestamp.hour), str(timestamp.minute), str(timestamp.second)
    os.chdir("%s" % ap)
    if "NULL" in prefix:
        os.system("mv dup_matrix.txt %s_dup_matrix.txt" % "".join(rename))
        os.system("mv names.txt %s_names.txt" % "".join(rename))
        os.system("mv duplicate_ids.txt %s_duplicate_ids.txt" % "".join(rename))
        os.system("mv paralog_ids.txt %s_paralog_ids.txt" % "".join(rename))
        os.system("mv bsr_matrix_values.txt %s_bsr_matrix.txt" % "".join(rename))
        if os.path.isfile("consensus.fasta"):
            os.system("mv consensus.fasta %s_consensus.fasta" % "".join(rename))
        if os.path.isfile("consensus.pep"):
            os.system("mv consensus.pep %s_consensus.pep" % "".join(rename))
    else:
        os.system("mv dup_matrix.txt %s_dup_matrix.txt" % prefix)
        os.system("mv names.txt %s_names.txt" % prefix)
        os.system("mv duplicate_ids.txt %s_duplicate_ids.txt" % prefix)
        os.system("mv paralog_ids.txt %s_paralog_ids.txt" % prefix)
        os.system("mv bsr_matrix_values.txt %s_bsr_matrix.txt" % prefix)
        if os.path.isfile("consensus.fasta"):
            os.system("mv consensus.fasta %s_consensus.fasta" % prefix)
        if os.path.isfile("consensus.pep"):
            os.system("mv consensus.pep %s_consensus.pep" % prefix)
    if "NULL" in prefix:
        outfile = open("%s_run_parameters.txt" % "".join(rename), "w")
    else:
        outfile = open("%s_run_parameters.txt" % prefix, "w")
    print >> outfile, "-d %s \\" % directory
    print >> outfile, "-i %s \\" % id
    print >> outfile, "-f %s \\" % filter
    print >> outfile, "-p %s \\" % processors
    print >> outfile, "-g %s \\" % genes
    print >> outfile, "-c %s \\" % cluster_method
    print >> outfile, "-b %s \\" % blast
    print >> outfile, "-l %s \\" % length
    print >> outfile, "-m %s \\" % max_plog
    print >> outfile, "-n %s \\" % min_hlog
    print >> outfile, "-t %s \\" % f_plog
    print >> outfile, "-k %s \\" % keep
    print >> outfile, "-s %s \\" % filter_peps
    print >> outfile, "-e %s \\" % filter_scaffolds
    print >> outfile, "-x %s \\" % prefix
    print >> outfile, "-z %s" % debug
    print >> outfile, "temp data stored here if kept: %s" % fastadir
    outfile.close()
    logging.logPrint("all Done")
    if "T" == keep:
        pass
    else:
        os.system("rm -rf %s" % fastadir)
    os.chdir("%s" % ap)

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-d", "--directory", dest="directory",
                      help="/path/to/fasta_directory [REQUIRED]",
                      type="string", action="callback", callback=test_dir)
    parser.add_option("-i", "--identity", dest="id", action="callback", callback=test_id,
                      help="clustering id threshold (0.0-1.0), defaults to 0.9",
                      type="float", default="0.9")
    parser.add_option("-f", "--filter", dest="filter", action="callback", callback=test_filter,
                      help="to use blast filtering or not, default is F or filter, change to T to turn off filtering",
                      default="F", type="string")
    parser.add_option("-p", "--parallel_workers", dest="processors",
                      help="How much work to do in parallel, defaults to 2",
                      default="2", type="int")
    parser.add_option("-g", "--genes", dest="genes", action="callback", callback=test_file,
                      help="predicted genes (nucleotide) to screen against genomes, will not use prodigal, must end in fasta (nt) or pep (aa)",
                      type="string",default="null")
    parser.add_option("-c", "--cluster_method", dest="cluster_method", action="callback", callback=test_cluster,
                      help="Clustering method to use: choose from usearch, vsearch, cd-hit",
                      type="string", default="NULL")
    parser.add_option("-b", "--blast", dest="blast", action="callback", callback=test_blast,
                      help="use tblastn, blastn, or blat (nucleotide search only), default is tblastn",
                      default="tblastn", type="string")
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
    parser.add_option("-s", "--filter_short_peps", dest="filter_peps", action="callback",
                      help="remove short peptides, smaller than 50AA?  Defaults to T",
                      default="T", callback=test_filter, type="string")
    parser.add_option("-e", "--filter_scaffolds", dest="filter_scaffolds", action="callback",
                      help="Filter any contig that contains an N? Defaults to F",
                      default="F", callback=test_filter, type="string")
    parser.add_option("-x", "--prefix", dest="prefix", action="store",
                      help="prefix for naming output files, defaults to time/date",
                      default="NULL", type="string")
    parser.add_option("-y", "--temp_dir", dest="temp_dir", action="store",
                      help="full path to desired temp directory location, defaults to python tempdir",
                      default="NULL", type="string")
    parser.add_option("-z", "--debug", dest="debug", action="callback",
                      help="turn debug on?  Defaults to F",
                      default="F", callback=test_filter, type="string")
    options, args = parser.parse_args()

    mandatories = ["directory"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.directory,options.id,options.filter,options.processors,options.genes,options.cluster_method,options.blast,
         options.length,options.max_plog,options.min_hlog,options.f_plog,options.keep,options.filter_peps,
         options.filter_scaffolds,options.prefix,options.temp_dir,options.debug)
