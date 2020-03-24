#!/usr/bin/env python

"""transfer annotation onto centroids"""

from __future__ import division
from __future__ import print_function
import optparse
import sys
import os
import subprocess
from collections import OrderedDict
from Bio import SeqIO
from ls_bsr.util import find_data_type

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s file cannot be opened' % option)
        sys.exit()

def blast_against_self(blast_type, query, database, output, processors):
    devnull = open('/dev/null', 'w')
    cmd = ["%s" % blast_type,
           "-query", query,
           "-db", database,
           "-num_threads", str(processors),
           "-evalue", "0.1",
           "-outfmt", "6",
           "-comp_based_stats", "F",
           "-seg", "no",
           "-out", output]
    subprocess.call(cmd, stdout=devnull, stderr=devnull)

def transfer_annotation(consensus, blast_in):
    blast_dict = {}
    outfile = open("genes.associated.fasta", "w")
    with open(blast_in) as infile:
        for line in infile:
            fields = line.split()
            blast_dict.update({fields[0]:fields[1]})
    with open(consensus) as in_fasta:
        for record in SeqIO.parse(in_fasta, "fasta"):
            if record.id in blast_dict:
                outfile.write(">"+str(blast_dict.get(record.id))+"\n")
                outfile.write(str(record.seq)+"\n")
            else:
                if int(len(record.seq))>10:
                    outfile.write(">"+record.id+"\n")
                    outfile.write(str(record.seq)+"\n")
    outfile.close()

def parse_blast_report(infile):
    """parse out only the name and bit score from the blast report"""
    outfile = open("query_blast.filtered", "w")
    with open(infile) as my_file:
        for line in my_file:
            try:
                fields = line.split("\t")
                outfile.write(fields[0]+"\t"+fields[1]+"\t"+fields[11])
            except:
                raise TypeError("malformed blast line found")
    outfile.close()

def get_unique_lines(infile):
    """only return the top hit for each query"""
    outfile = open("query.filtered.unique", "w")
    d = {}
    with open(infile) as my_file:
        for line in my_file:
            unique = line.split("\t",1)[0]
            if unique not in d:
                d[unique] = 1
                outfile.write(line,)
    outfile.close()

def get_cluster_ids(in_fasta):
    clusters = []
    with open(in_fasta) as infile:
        for record in SeqIO.parse(infile, "fasta"):
            clusters.append(record.id)
    nr = list(OrderedDict.fromkeys(clusters))
    if len(clusters) == len(nr):
        return clusters
    else:
        print("Problem with gene list.  Are there duplicate headers in your file?")
        sys.exit()

def parse_self_blast(lines):
    my_dict={}
    for line in lines:
        try:
            fields=line.split()
            str1=fields[0]
            str2=fields[11]
            my_dict.update({str1:str2})
        except:
            raise TypeError("blast file is malformed")
    return my_dict

def process_consensus(in_fasta,new_dict,out_fasta_prefix):
    my_lists = []
    outfile = open("%s.annotated.fasta" % out_fasta_prefix, "w")
    with open(in_fasta) as my_fasta:
        for record in SeqIO.parse(my_fasta, "fasta"):
            if record.id not in new_dict:
                outfile.write(">"+str(record.id)+"\n")
                outfile.write(str(record.seq)+"\n")
            else:
                new_id = record.id.replace(record.id, new_dict.get(record.id))
                outfile.write(">"+str(new_id)+"::"+record.id+"\n")
                outfile.write(str(record.seq)+"\n")
    outfile.close()


def update_dict(ref_scores,query_file,all_clusters,threshold):
    new_dict = {}
    with open(query_file) as my_file:
        for line in my_file:
            newline = line.strip()
            fields = newline.split()
            """Orders by cluster"""
            for cluster in all_clusters:
                if cluster == fields[1]:
                    try:
                        if (float(fields[2])/float(ref_scores.get(fields[1]))*100)>int(threshold):
                            new_dict.update({fields[0]:fields[1]})
                    except:
                        print("couldn't process", fields[2], ref_scores.get(fields[0]), fields[1], fields[0])
    """Returns centroid:associated_gene"""
    return new_dict,len(new_dict)

def main(aligner,peptides,consensus,processors,threshold,out_fasta_prefix):
    devnull = open("/dev/null", "w")
    """These are your reference peptides"""
    pep_path=os.path.abspath("%s" % peptides)
    consensus_path=os.path.abspath("%s" % consensus)
    data_type = find_data_type(consensus_path)
    #data_type is either aa pr mt
    if consensus_path.endswith(".pep"):
        if data_type == "aa":
            pass
        else:
            print("Make sure that nucleotide data ends in 'fasta'")
            sys.exit()
        if aligner == "blastp":
            blast_type = "blastp"
            ab = subprocess.call(['which', 'blastp'])
            if ab == 0:
                pass
            else:
                print("blastp must be in your path to use peptides")
                sys.exit()
        else:
            blast_type = "diamond"
            ab = subprocess.call(['which', 'diamond'])
            if ab == 0:
                pass
            else:
                print("diamond is not in your path but needs to be")
                sys.exit()
    elif consensus_path.endswith(".fasta"):
        """Need to check to make sure that this is a nucleotide file"""
        if data_type == "nt":
            pass
        else:
            print("Make sure that peptide data ends with 'pep'")
            sys.exit()
        blast_type = "blastx"
        ab = subprocess.call(['which', 'blastx'])
        if ab == 0:
            pass
        else:
            print("blastx must be in your path to use nucleotides")
            sys.exit()
    """"This causes problems"""
    os.system("cp %s query.peptides.xyx" % pep_path)
    if blast_type == "diamond":
        subprocess.check_call("diamond makedb --in query.peptides.xyx -d DB > /dev/null 2>&1", stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
    else:
        subprocess.check_call("makeblastdb -in query.peptides.xyx -dbtype prot > /dev/null 2>&1", stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
    if blast_type == "blastp":
        blast_against_self("blastp", "query.peptides.xyx", "query.peptides.xyx", "xyx.blast.out", processors)
    elif blast_type == "diamond":
        subprocess.check_call("diamond blastp -p %s -d DB -f 6 -q query.peptides.xyx -o xyx.blast.out > /dev/null 2>&1" % processors, stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
    else:
        blast_against_self("blastx", "query.peptides.xyx", "query.peptides.xyx", "xyx.blast.out", processors)
    os.system("sort -u -k 1,1 xyx.blast.out > xyx.blast.unique.xyx")
    ref_scores=parse_self_blast(open("xyx.blast.unique.xyx"))
    ref_score_file = open("ref.scores", "w")
    for k,v in ref_scores.items():
        ref_score_file.write(str(k)+"\t"+str(v)+"\n")
    ref_score_file.close()
    os.system("rm xyx.blast.out xyx.blast.unique.xyx")
    if blast_type == "blastp":
        blast_against_self("blastp", consensus_path, "query.peptides.xyx", "xyx.blast.out", processors)
    elif blast_type == "blastx":
        blast_against_self("blastx", consensus_path, "query.peptides.xyx", "xyx.blast.out", processors)
    elif blast_type == "diamond":
        subprocess.check_call("diamond blastp -p %s -d DB -f 6 -q %s -o xyx.blast.out > /dev/null 2>&1" % (processors,consensus_path), stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
    parse_blast_report("xyx.blast.out")
    get_unique_lines("query_blast.filtered")
    """These are reference clusters"""
    clusters = get_cluster_ids(pep_path)
    """these are query clusters"""
    ref_clusters = get_cluster_ids(consensus_path)
    """new dict should only include those pairs were annotation was transferred"""
    new_dict, dict_len = update_dict(ref_scores, "query.filtered.unique", clusters, threshold)
    process_consensus(consensus_path,new_dict,out_fasta_prefix)
    missing = []
    for ref_cluster in ref_clusters:
        if ref_cluster not in new_dict:
            missing.append(ref_cluster)
    print("Total number of query peptides = %s" % len(ref_clusters))
    print("Number of query peptides with transferred annotation = %s" % dict_len)
    if len(missing)>0:
        outfile = open("missing_ids.txt", "w")
        print("Number of query peptides with no transferred annotation = %s" % len(missing))
        for loci in missing:
            outfile.write(str(loci)+"\n")
        outfile.close()
    os.system("rm query.peptides.xyx* xyx.blast.out query_blast.filtered query.filtered.unique")

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-a", "--aligner", dest="aligner",
                      help="aligner to use; choose from blastp or diamond (default)",
                      type="string", action="store", default="diamond")
    parser.add_option("-p", "--peptides", dest="peptides",
                      help="/path/to/annotated peptides (reference) [REQUIRED]",
                      type="string", action="callback", callback=test_file)
    parser.add_option("-c", "--consensus", dest="consensus", action="callback", callback=test_file,
                      help="/path/to/consensus file, can be nucleotide or peptide [REQUIRED]",
                      type="string")
    parser.add_option("-r", "--processors", dest="processors", action="store",
                      help="number of processors to use with BLAST, defaults to 2",
                      type="int", default="2")
    parser.add_option("-t", "--threshold", dest="threshold",
                      help="[integer] lower percentage threshold for assigning annotation, defaults to 80[%]",
                      type="int", action="store", default="80")
    parser.add_option("-o", "--out_fasta_prefix", dest="out_fasta_prefix",
                      help="Naming scheme for output files [REQUIRED]",
                      type="string", action="store")
    options, args = parser.parse_args()

    mandatories = ["peptides", "consensus", "out_fasta_prefix"]
    for m in mandatories:
        if not getattr(options, m, None):
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)

    main(options.aligner,options.peptides,options.consensus,options.processors,options.threshold,options.out_fasta_prefix)
