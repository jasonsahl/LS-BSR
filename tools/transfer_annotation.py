#!/usr/bin/env python

"""transfer annotation onto centroids"""

import optparse
import sys
import os
import subprocess
from Bio import SeqIO

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def blast_against_self_tblastn(blast_type, query, database, output, processors):
    devnull = open('/dev/null', 'w')
    cmd = ["%s" % blast_type,
           "-query", genes_pep,
           "-db", genes_nt,
           "-num_threads", str(processors),
           "-evalue", "0.1",
           "-outfmt", "6",
           "-out", output]
    subprocess.call(cmd, stdout=devnull, stderr=devnull)


def transfer_annotation(consensus, blast_in):
    blast_dict = {}
    outfile = open("genes.associated.fasta", "w")
    for line in open(blast_in, "U"):
        fields = line.split()
        blast_dict.update({fields[0]:fields[1]})
    for record in SeqIO.parse(consensus, "fasta"):
        if record.id in blast_dict:
            print >> outfile,">"+str(blast_dict.get(record.id))
            print >> outfile,record.seq
        else:
            if int(len(record.seq))>10:
                print >> outfile,">"+record.id
                print >> outfile,record.seq
    outfile.close()

def parse_blast_report(infile):
    """parse out only the name and bit score from the blast report"""
    outfile = open("query_blast.filtered", "w")
    for line in open(infile, "rU"):
        try:
            fields = line.split("\t")
            print >> outfile, fields[0]+"\t"+fields[1]+"\t"+fields[11],
        except:
            raise TypeError("malformed blast line found")
    outfile.close()

def get_unique_lines(infile):
    """only return the top hit for each query"""
    outfile = open("query.filtered.unique", "w")
    d = {}
    input = file(infile)
    for line in input:
        unique = line.split("\t",1)[0]
        if unique not in d:
            d[unique] = 1
            print >> outfile,line,
    outfile.close()

def get_cluster_ids(in_fasta):
    clusters = []
    infile = open(in_fasta, "U")
    for record in SeqIO.parse(infile, "fasta"):
        clusters.append(record.id)
    nr = list(OrderedDict.fromkeys(clusters))
    if len(clusters) == len(nr):
        return clusters
    else:
        print "Problem with gene list.  Are there duplicate headers in your file?"
        sys.exit()

    
def main(peptides,consensus,processors,threshold):
    devnull = open("/dev/null", "w")
    pep_path=os.path.abspath("%s" % peptides)
    consensus_path=os.path.abspath("%s" % consensus)
    if consensus_path.endswith(".pep"):
        ab = subprocess.call(['which', 'blastp'])
        if ab == 0:
            pass
        else:
            print "blastp must be in your path to use peptides"
            sys.exit()
    elif consensus_path.endswith(".fasta"):
        ab = subprocess.call(['which', 'blastx'])
        if ab == 0:
            pass
        else:
            print "blastx must be in your path to use nucleotides"
            sys.exit()
    """"removes empty white space from your input file"""
    os.system("sed 's/ /_/g' %s > query.peptides.xyx" % pep_path)
    subprocess.check_call("makeblastdb -in query.peptides.xyx -dbtype prot", shell=True)
    if consensus_path.endswith(".pep"):
        blast_against_self("blastp", "query.peptides.xyx", "query.peptides.xyx", "xyx.blast.out", processors)
    elif consensus_path.endswith(".fasta"):
        blast_against_self("blastx", "query.peptides.xyx", "query.peptides.xyx", "xyx.blast.out", processors)
    else:
        print "input genes file is of incorrect format, choose from fasta or pep"
        sys.exit()
    os.system("sort -u -k 1,1 xyx.blast.out.xyx > xyx.blast.unique.xyx")
    ref_scores=parse_self_blast(open("xyx.blast.unique.xyx", "U"))
    os.system("rm xyx.blast.out xyx.blast.unique.xyx")
    if consensus_path.endswith(".pep"):
        blast_against_self("blastp", consensus_path, "query.peptides.xyx", "xyx.blast.out", processors)
    elif consensus_path.endswith(".fasta"):
        blast_against_self("blastx", consensus_path, "query.peptides.xyx", "xyx.blast.out", processors)
    parse_blast_report("xyx.blast.out")
    get_unique_lines("query_blast.filtered")
    clusters = get_cluster_ids(pep_path)
    new_dict = update_dict(ref_scores, "query.filtered.unique", clusters, threshold)
    process_consensus(consensus,new_dict)
    #transfer_annotation(consensus, "xyx.blast.unique.xyx")
    #os.system("rm query.peptides.xyx")

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-p", "--peptides", dest="peptides",
                      help="/path/to/annotated peptides [REQUIRED]",
                      type="string", action="callback", callback=test_file)
    parser.add_option("-c", "--consensus", dest="consensus", action="callback", callback=test_file,
                      help="/path/to/consensus file, can be nucleotide or peptide [REQUIRED]",
                      type="string")
    parser.add_option("-r", "--processors", dest="processors", action="store",
                      help="number of processors to use with BLAST, defaults to 2",
                      type="int", default="2")
    parser.add_option("-t", "--threshold", dest="threshold",
                      help="[integer] lower BSR threshold for assigning annotation, defaults to 80[%]",
                      type="int", action="store", default="80"
    options, args = parser.parse_args()
    
    mandatories = ["peptides", "consensus"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.peptides,options.consensus,options.processors,options.threshold)
    
