#!/usr/bin/env python

"""Calculate the BSR value for all predicted ORFs
in a set of genomes in fasta format.  V3 - replaced
transeq with BioPython.  V4 - changed to true BSR
test.  V5 - fixed bug in how BSR was calculated.
V6 - changed gene caller from Glimmer to Prodigal"""

import Bio
from Bio import SeqIO
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
import sys
import os
import glob
import optparse
import subprocess
import shlex
from subprocess import call
from Bio.SeqRecord import SeqRecord
from igs.utils import functional as func
from igs.utils import logging
from igs.threading import functional as p_func
import errno
import unittest
from Bio.Seq import Seq
import threading

rec=1

def autoIncrement(): 
    global rec 
    pStart = 1  
    pInterval = 1 
    if (rec == 0):  
        rec = pStart  
    else:  
        rec += pInterval  
        return rec

def get_seq_name(in_fasta):
    """used for renaming the sequences"""
    return os.path.basename(in_fasta)

def predict_genes(directory, processors):
    """simple gene prediction using Prodigal in order
    to find coding regions from a genome sequence"""    
    os.chdir("%s/joined" % directory)
    curr_dir=os.getcwd()
    files = os.listdir(curr_dir)
    files_and_temp_names = [(str(idx), os.path.join(curr_dir, f))
                            for idx, f in enumerate(files)]
                         
    def _perform_workflow(data):
        tn, f = data
        subprocess.check_call("prodigal -i %s -d %s_genes.seqs -a %s_genes.pep > /dev/null 2>&1" % (f, reduced, reduced), shell=True)
        
    results = set(p_func.pmap(_perform_workflow,
                              files_and_temp_names,
                              num_workers=processors))

def uclust_sort():
    """sort with Usearch. Updated to V6"""
    cmd = ["usearch6", 
           "-sortbylength", "all_gene_seqs.out",
           "-output", "all_sorted.txt"]
    subprocess.check_call(cmd)

def uclust_cluster(id):
    """cluster with Uclust.  Updated to V6"""
    cmd = ["usearch6",
           "-cluster_fast", "all_sorted.txt",
           "-id", str(id),
           "-uc", "results.uc",
           "-centroids", "consensus.fasta"]
    subprocess.check_call(cmd)

def translate_consensus(consensus):
    """translate nucleotide into peptide with BioPython"""
    infile = open(consensus, "rU")
    #output_handle = open("tmp.pep", "w")
    output = [ ]
    for record in SeqIO.parse(infile, "fasta"):
        #print >> output_handle, ">"+record.id
        #print >> output_handle, record.seq.translate(to_stop=True)
        output.append(record.seq.translate(to_stop=True))
    for record in output: return str(record)

def filter_seqs(input_pep):
    """filter out short sequences from a multifasta.
    Will hopefully speed up the process without losing
    important information"""
    long_sequences = [ ]
    infile = open(input_pep, "rU")
    #outfile = open("consensus.pep", "w")
    output = [ ]
    for record in SeqIO.parse(infile, "fasta"):
        if len(record.seq) > int(50):
            long_sequences.append(record)
            output.append(len(record.seq))
            #SeqIO.write(long_sequences, outfile, "fasta")
            #outfile.close()
    return output

def blast_against_each_genome(directory, processors, filter):
    """BLAST all peptides against each genome"""
    curr_dir=os.getcwd()
    files = os.listdir(curr_dir)
    files_and_temp_names = [(str(idx), os.path.join(curr_dir, f))
                            for idx, f in enumerate(files)]
    def _perform_workflow(data):
	tn, f = data
        if ".fasta.new" in f:
            subprocess.check_call("formatdb -i %s -p F > /dev/null 2>&1" % f, shell=True)
        if ".fasta.new" in f:
            cmd = ["blastall",
                   "-p", "tblastn",
                   "-i", "consensus.pep",
                   "-d", f,
                   "-a", str(processors),
                   "-e", "0.1",
                   "-m", "8",
                   "-F", str(filter),
                   "-o", "%s_blast.out" % f]
            subprocess.check_call(cmd)
            
    results = set(p_func.pmap(_perform_workflow,
                              files_and_temp_names,
                              num_workers=processors))

def parse_blast_report(directory):
    """parse out only the name and bit score from the blast report"""
    curr_dir=os.getcwd()
    for infile in glob.glob(os.path.join(curr_dir, "*_blast.out")):
        names = get_seq_name(infile)
        ref = open(infile, "rU")
        data = ref.readlines()
        #outfile = open("%s.filtered" % names, "w")
        outdata = [ ]
        for line in data:
            fields = line.split()
            #print >> outfile, fields[0]+"\t"+fields[11],
            outdata.append(fields[0])
            outdata.append(fields[11])
        return outdata
            
def get_unique_lines(directory):
    """only return the top hit for each query"""
    curr_dir=os.getcwd()
    for infile in glob.glob(os.path.join(curr_dir, "*.filtered")):
        names = get_seq_name(infile)
        #outfile = open("%s.filtered.unique" % names, "w")
        d = {}
        input = file(infile)
        outdata = [ ]
        for line in input:
            unique = line.split("\t",1)[0]
            if unique not in d:
                d[unique] = 1
                #print >> outfile,line,
                outdata.append(line)
        return outdata
    
def make_table(directory, processors):
    """make the BSR matrix table"""
    clusters=[ ]
    curr_dir=os.getcwd()
    """I only use this loop to grab names...combine with next loop?
       I need the nr values before the next loop"""
    for infile in glob.glob(os.path.join(curr_dir, "*.filtered.unique")):
        file=open(infile, "rU")
        for line in file:
		fields=line.split()
                if fields[0] not in clusters:
                    clusters.append(fields[0])
    """de-replicate the clusters"""
    nr=[x for i, x in enumerate(clusters) if x not in clusters[i+1:]]
    names = [ ]
    outdata = [ ]
    files = glob.glob(os.path.join(curr_dir, "*.filtered.unique"))
    files_and_temp_names = [(str(idx), os.path.join(curr_dir, f))
                            for idx, f in enumerate(files)]
    lock = threading.Lock()
    def _perform_workflow(data):
        lock.acquire()
        tn, f = data
        """get the name of each of the files to be iterated"""
        name=[ ]
        out=get_seq_name(f)
        name.append(out)
        reduced=[ ]
        """remove the junk at the end of the file"""
        for x in name:reduced.append(x.replace('.filtered.unique',''))
        names.append(reduced)
        dict={}
        file=open(f, "rU")
        #tmpfile=open("tmp.txt", "w")
        """make a dictionary of all clusters and values"""
        for line in file:
            fields=line.split()
            dict.update({fields[0]:fields[1]})
        cluster_names={}
        """add in values, including any potentially missing ones"""
        for k,v in dict.iteritems():
            if k in nr: cluster_names.update({k:v})
        for x in nr:
            if x not in dict.keys():cluster_names.update({x:0})
        """need to write a blank space"""
                #for x in reduced: open("%s.tmp.matrix" % x, 'a').write('%s\n' % x)
        """sort keys to get the same order between samples"""
        for key in sorted(cluster_names.iterkeys()):
            for x in reduced:
                outdata.append(cluster_names[key])
                #open("%s.tmp.matrix" % x, 'a').write("%s\n" % cluster_names[key])
        lock.release()
        #print cluster_names
        #print reduced
    results = set(p_func.pmap(_perform_workflow,
                              files_and_temp_names,
                              num_workers=processors))
        #names_out = open("names.txt", "w")
        #for x in names: print >> names_out, "".join(x)
        #nr_sorted=sorted(nr)
    #open("ref.list", "a").write("\n")
    #for x in nr_sorted:
    #open("ref.list", "a").write("%s\n" % x)
    return outdata


def divide_values(file, ref_scores):
    """divide each BSR value in a row by that row's maximum value"""
    infile = open(file, "rU")
    firstLine = infile.readline()
    FL_F=firstLine.split()
    #outfile = open("BSR_matrix_values.txt", "a")
    #print >> outfile, '\t'.join([str(item) for item in FL_F])
    outdata=[ ]
    for line in infile:
        fields=line.split()
        all_fields=list(fields)
        fields=map(float, fields[1:])
        #print fields
        values= [ ]
        """v7 - added in exception to address bug"""
	for x in fields:
            try:
                values.append(float(x)/float(ref_scores.get(all_fields[0])))
            except:
                values.append(float(x)/float("1000"))
        sort_values=['%.2f' % elem for elem in values]
        #print >> outfile, '\t'.join([str(item) for item in sort_values])
        outdata.append(values)
    return outdata

def blast_against_self(genes_nt, genes_pep, output, filter):
    cmd = ["blastall",
           "-p", "tblastn",
           "-i", genes_pep,
           "-d", genes_nt,
           "-a", "2",
           "-e", "0.1",
           "-m", "8",
           "-F", str(filter),
           "-o", output]
    subprocess.check_call(cmd)

def parse_self_blast(blast_out):
    my_dict={}
    for line in open(blast_out):
        fields=line.split()
        str1=fields[0]
        str2=fields[11]
        my_dict.update({str1:str2})
    return my_dict

def translate_genes(genes, directory):
    """translate nucleotide into peptide with BioPython"""
    #os.chdir("%s/joined" % directory)
    infile = open(genes, "rU")
    output = [ ]
    #output_handle = open("genes.pep", "w")
    for record in SeqIO.parse(infile, "fasta"):
        if len(record.seq.translate(to_stop=True, table=11))>=30:
            output.append(record.seq.translate(to_stop=True, table=11))
            #print >> output_handle, ">"+record.id
            #print >> output_handle, record.seq.translate(to_stop=True, table=11)
    for record in output: return str(record)
    infile.close()
    #output_handle.close()

def rename_fasta_header(fasta_in):
    """this is used for parsing the mugsy output"""
    #handle = open(fasta_out, "w")
    outdata = [ ]
    for record in SeqIO.parse(open(fasta_in), "fasta"):
        outdata.append(">"+"centroid"+"_"+str(autoIncrement()))
        #print >> handle, ">"+"centroid"+"_"+str(autoIncrement())
        #print >> handle, record.seq
    return outdata
    handle.close()

class Test1(unittest.TestCase):
    def test(self):
        self.assertEqual(get_seq_name("/path/to/test.fasta"), "test.fasta")
class Test2(unittest.TestCase):
    def test(self):
        self.assertEqual(translate_consensus("test_consensus.fasta"), 'MTSFP')
class Test3(unittest.TestCase):
    def test(self):
        self.assertEqual(filter_seqs("test_length.pep"), [62])
class Test4(unittest.TestCase):
    def test(self):
        self.assertEqual(parse_self_blast("blast.test"), {'Cluster2': '60.6', 'Cluster0': '30.2', 'Cluster1': '40.5'})
class Test5(unittest.TestCase):
    def test(self):
        self.assertEqual(parse_blast_report("."), ['Cluster0', '30.2'])
class Test6(unittest.TestCase):
    def test(self):
        self.assertEqual(get_unique_lines("."), ['Cluster0\t30.2'])
class Test7(unittest.TestCase):
    def test(self):
        self.assertEqual(make_table(".", 2), ['30.2', '40.5', '60.6', '15.2', 0, '30.6'])
class Test8(unittest.TestCase):
    def test(self):
        self.assertEqual(divide_values( "test_bsr.matrix", {'Cluster2': '60.6', 'Cluster0': '30.2', 'Cluster1': '40.5'}),
                     [[1.0, 0.5033112582781457], [1.0, 0.0], [1.0, 0.504950495049505]])
class Test9(unittest.TestCase):
    def test(self):
        self.assertEqual(translate_genes("genes_test.fasta", "."), 'MNPHLTEHPPVGDIDALLQDTWLQVISLRQGVTCAEGEGQAFWQRCVADIERVHQALKDAGHSEQSCQHIRYAQCALLDETVKGRGVQDDAYFVWCHSPLQAHFFNTLDAGSQLYERMRAVLREPAPDRAVLTCFHRVLMLGFLGGYASPAASEREQLIDQLSVQVPAFSVAPSRGILASAASRNRLGIWLRYWPVRLGLAALMVALLWWGLDHWLSGLLATLLPEPV')
class Test10(unittest.TestCase):
    def test(self):
        self.assertEqual(rename_fasta_header("genes_test.fasta"), ['>centroid_2', '>centroid_3'])
class Test11(unittest.TestCase):
    def test(self):
        self.assertEqual(autoIncrement(), 4)


if __name__ == "__main__":
    unittest.main()
   

