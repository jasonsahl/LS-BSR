#!/usr/bin/env python

"""annotate a BSR matrix,
based on a GenBank file"""


import optparse
import sys
import os
from Bio import SeqIO
import subprocess
from collections import OrderedDict
from operator import itemgetter

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def blat_against_self(query,reference,output):
    subprocess.check_call("blat -out=blast8 -minIdentity=75 %s %s %s > /dev/null 2>&1" % (reference,query,output), shell=True)

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

def update_dict(ref_scores, query_file, all_clusters, threshold):
    new_dict = {}
    for line in open(query_file, "U"):
        newline = line.strip()
        fields = newline.split()
        hits = []
        for cluster in all_clusters:
            if cluster == fields[1]:
                if (float(fields[2])/float(ref_scores.get(fields[0]))*100)>int(threshold):
                    new_dict.update({fields[1]:fields[0]})
                    hits.append("1")
        if len(hits) == 0:
            new_dict.update({fields[0]:fields[0]})
    for cluster in all_clusters:
        if cluster in new_dict:
            pass
        else:
            new_dict.update({cluster:cluster})
    return new_dict

def process_bsr_matrix(matrix,new_dict):
    my_lists = []
    outfile = open("bsr_matrix_annotated.txt", "w")
    for line in open(matrix, "U"):
        if not line.startswith("centroid"):
            refline = line.strip()
            ref_fields = refline.split()
            ref_fields.insert(0,"")
            my_lists.append(ref_fields)
        else:
            newline = line.strip()
            fields = newline.split()
            fields[0] = new_dict.get(fields[0])
            my_lists.append(fields)
    sorted_lists = sorted(my_lists, key=itemgetter(0))
    for alist in sorted_lists:
        print >> outfile, "\t".join(alist)
    outfile.close()

def process_consensus(consensus,new_dict):
    outfile = open("consensus_annotated.fasta", "w")
    for record in SeqIO.parse(open(consensus, "U"), "fasta"):
        changed_id = []
        for k,v in new_dict.iteritems():
            if k == record.id:
                new_id = record.id.replace("%s" % k, "%s" % v)
                changed_id.append(new_id)
        if len(changed_id)>0:
            print >> outfile, ">"+"".join(changed_id)+"\n"+str(record.seq)
        else:
            print >> outfile, ">"+str(record.id)+"\n"+str(record.seq)
    outfile.close()

def get_seq_name(in_fasta):
    """used for renaming the sequences"""
    return os.path.basename(in_fasta)

def main(bsr_matrix, consensus, locus_tags, threshold):
    ac = subprocess.call(['which', 'blat'])
    if ac == 0:
        pass
    else:
        print "You have requested blat, but it is not in your PATH"
        sys.exit()
    name = get_seq_name(locus_tags)
    locus_path = os.path.abspath("%s" % locus_tags)
    blat_against_self("%s" % locus_path, "%s" % locus_path, "tmp_blast.out")
    subprocess.check_call("sort -u -k 1,1 tmp_blast.out > self_blast.out", shell=True)
    ref_scores=parse_self_blast(open("self_blast.out", "U"))
    blat_against_self("%s" % locus_path, "%s" % consensus, "query_blast.out")
    parse_blast_report("query_blast.out")
    get_unique_lines("query_blast.filtered")
    clusters = get_cluster_ids(consensus)
    new_dict = update_dict(ref_scores, "query.filtered.unique", clusters, threshold)
    process_bsr_matrix(bsr_matrix,new_dict)
    process_consensus(consensus,new_dict)
    os.system("rm tmp_blast.out self_blast.out query_blast.out query_blast.filtered query.filtered.unique")

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-b", "--bsr_matrix", dest="bsr_matrix",
                      help="/path/to/bsr_matrix [REQUIRED]",
                      type="string", action="callback", callback=test_file)
    parser.add_option("-c", "--consensus", dest="consensus",
                      help="/path/to/consensus.fasta file (nucleotides) [REQUIRED]",
                      type="string", action="callback", callback=test_file)
    parser.add_option("-l", "--locus_tags", dest="locus_tags",
                      help="/path/to/locus_tags file (nucleotides) [REQUIRED]",
                      type="string", action="callback", callback=test_file)
    parser.add_option("-t", "--threshold", dest="threshold",
                      help="[integer] lower BSR threshold for assigning annotation, defaults to 80[%]",
                      type="int", action="store", default="80")
    options, args = parser.parse_args()

    mandatories = ["bsr_matrix", "consensus", "locus_tags"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.bsr_matrix,options.consensus,options.locus_tags,options.threshold)
