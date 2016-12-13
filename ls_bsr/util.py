#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
import sys
import os
import glob
import optparse
import subprocess
from subprocess import Popen
import shlex
from subprocess import call
import random
import collections
try:
    from Bio.SeqRecord import SeqRecord
    import Bio
    from Bio import SeqIO
    from Bio import Phylo
except:
    print("BioPython is not in your PATH, but needs to be")
    sys.exit()
try:
    from igs.utils import functional as func
    from igs.utils import logging
    from igs.threading import functional as p_func
    """test code"""
except:
    print("Your environment is not set correctly.  Please add LS-BSR to your PYTHONPATH and try again")
    sys.exit()
import errno
import threading
import types
from collections import deque,OrderedDict
import collections

def get_cluster_ids(in_fasta):
    clusters = []
    infile = open(in_fasta, "U")
    for record in SeqIO.parse(infile, "fasta"):
        clusters.append(record.id)
    nr = list(OrderedDict.fromkeys(clusters))
    if len(clusters) == len(nr):
        return clusters
    else:
        print("Problem with gene list.  Are there duplicate headers in your file?")
        sys.exit()

def divide_values(file, ref_scores):
    """divide each BSR value in a row by that row's maximum value"""
    errors = []
    infile = open(file, "rU")
    firstLine = infile.readline()
    FL_F=firstLine.split()
    outfile = open("BSR_matrix_values.txt", "w")
    outfile.write('\t'.join([str(item) for item in FL_F])+"\n")
    outdata=[]
    for line in infile:
        fields=line.split()
        all_fields=list(fields)
        try:
            fields=map(float, all_fields[1:])
        except:
            raise TypeError("abnormal number of fields observed")
        values= [ ]
        for x in fields:
            try:
                values.append(float(x)/float(ref_scores.get(all_fields[0])))
            except:
                """if a mismatch error in names encountered, change values to 0"""
                errors.append(all_fields[0])
                values.append(float("0"))
        sort_values=['%.2f' % elem for elem in values]
        outfile.write('\t'.join([str(item) for item in sort_values])+"\n")
        outdata.append(values)
    if len(errors)>0:
        nr=[x for i, x in enumerate(errors) if x not in errors[i+1:]]
        logging.logPrint("The following genes had no hits in datasets or are too short, values changed to 0, check names and output: %s" % "\n".join(nr))
    outfile.close()
    return outdata

def predict_genes(fastadir, processors):
    """simple gene prediction using Prodigal in order
    to find coding regions from a genome sequence"""
    os.chdir("%s" % fastadir)
    files = os.listdir(fastadir)
    files_and_temp_names = [(str(idx), os.path.join(fastadir, f))
                            for idx, f in enumerate(files)]
    def _perform_workflow(data):
        tn, f = data
        subprocess.check_call("prodigal -i %s -d %s_genes.seqs -a %s_genes.pep > /dev/null 2>&1" % (f, f, f), shell=True)
    results = set(p_func.pmap(_perform_workflow,
                              files_and_temp_names,
                              num_workers=processors))

def predict_genes_dev(fastadir, processors):
    """simple gene prediction using Prodigal in order
    to find coding regions from a genome sequence"""
    os.chdir("%s" % fastadir)
    files = os.listdir(fastadir)
    files_and_temp_names = [(str(idx), os.path.join(fastadir, f))
                            for idx, f in enumerate(files)]
    def _perform_workflow(data):
        tn, f = data
        subprocess.check_call("prodigal -i %s -d %s_genes.seqs -a %s_genes.pep > /dev/null 2>&1" % (f, f, f), shell=True)
    results = mp_shell(_perform_workflow, files_and_temp_names, processors)


def rename_fasta_header(fasta_in, fasta_out):
    """this is used for renaming the output,
    in the off chance that there are duplicate
    names for separate peptides"""
    rec=1
    handle = open(fasta_out, "w")
    outdata = [ ]
    for record in SeqIO.parse(open(fasta_in), "fasta"):
        try:
            outdata.append(">"+"centroid"+"_"+record.id)
            handle.write(">centroid_"+str(autoIncrement())+"\n")
            handle.write(str(record.seq)+"\n")
        except:
            raise TypeError("problem with input sequence encountered")
    handle.close()
    return outdata

def uclust_cluster(id):
    devnull = open("/dev/null", "w")
    """cluster with Uclust.  Updated to V6"""
    cmd = ["usearch",
           "-cluster_fast", "all_sorted.txt",
           "-id", str(id),
           "-uc", "results.uc",
           "-centroids", "consensus.fasta"]
    subprocess.call(cmd, stderr=devnull, stdout=devnull)

def blast_against_each_genome_tblastn(dir_path, processors, peptides, filter):
    """BLAST all peptides against each genome"""
    curr_dir=os.getcwd()
    files = os.listdir(curr_dir)
    devnull = open("/dev/null", "w")
    if "T" in filter:
        my_seg = "yes"
    else:
        my_seg = "no"
    files_and_temp_names = [(str(idx), os.path.join(curr_dir, f))
                            for idx, f in enumerate(files)]
    def _perform_workflow(data):
        tn, f = data
        if ".fasta.new" in f:
            try:
                subprocess.check_call("makeblastdb -in %s -dbtype nucl > /dev/null 2>&1" % f, shell=True)
            except:
                print("problem found in formatting genome %s" % f)
        if ".fasta.new" in f:
            try:
                devnull = open('/dev/null', 'w')
                cmd = ["tblastn",
                       "-query", peptides,
                       "-db", f,
                       "-seg", my_seg,
                       "-comp_based_stats", "F",
                       "-num_threads", "1",
                       "-evalue", "0.1",
                       "-outfmt", "6",
                       "-out", "%s_blast.out" % f]
                subprocess.call(cmd, stdout=devnull, stderr=devnull)
            except:
                print("genomes %s cannot be used" % f)

    results = set(p_func.pmap(_perform_workflow,
                              files_and_temp_names,
                              num_workers=processors))

def blast_against_each_genome_blastn(dir_path, processors, filter, peptides):
    """BLAST all peptides against each genome"""
    if "F" in filter:
        my_seg = "yes"
    else:
        my_seg = "no"
    curr_dir=os.getcwd()
    files = os.listdir(curr_dir)
    files_and_temp_names = [(str(idx), os.path.join(curr_dir, f))
                            for idx, f in enumerate(files)]
    def _perform_workflow(data):
        tn, f = data
        if ".fasta.new" in f:
            try:
                subprocess.check_call("makeblastdb -in %s -dbtype nucl > /dev/null 2>&1" % f, shell=True)
            except:
                print("problem found in formatting genome %s" % f)
        if ".fasta.new" in f:
            devnull = open('/dev/null', 'w')
            try:
                cmd = ["blastn",
                       "-task", "blastn",
                       "-query", peptides,
                       "-db", f,
                       "-dust", str(my_seg),
                       "-num_threads", "1",
                       "-evalue", "0.1",
                       "-outfmt", "6",
                       "-out", "%s_blast.out" % f]
                subprocess.call(cmd, stdout=devnull, stderr=devnull)
            except:
                print("The genome file %s was not processed" % f)

    results = set(p_func.pmap(_perform_workflow,
                              files_and_temp_names,
                              num_workers=processors))

def get_seq_name(in_fasta):
    """used for renaming the sequences"""
    return os.path.basename(in_fasta)

def filter_seqs(input_pep):
    """filter out short sequences from a multifasta.
    Will hopefully speed up the process without losing
    important information"""
    long_sequences = [ ]
    infile = open(input_pep, "rU")
    outfile = open("consensus.pep", "w")
    outdata = [ ]
    for record in SeqIO.parse(infile, "fasta"):
        if len(record.seq) >= int(50):
            long_sequences.append(record)
            outdata.append(len(record.seq))
    SeqIO.write(long_sequences, outfile, "fasta")
    outfile.close()
    return outdata

def parse_blast_report(test):
    """parse out only the name and bit score from the blast report"""
    curr_dir=os.getcwd()
    outdata = [ ]
    for infile in glob.glob(os.path.join(curr_dir, "*_blast.out")):
        names = get_seq_name(infile)
        outfile = open("%s.filtered" % names, "w")
        for line in open(infile, "rU"):
            try:
                fields = line.split("\t")
                outfile.write(fields[0]+"\t"+fields[11],)
                if "true" in test:
                    outdata.append(fields[0])
                    outdata.append(fields[11])
                else:
                    pass
            except:
                raise TypeError("malformed blast line found")
        outfile.close()
    if "true" in test:
        return outdata
    else:
        pass

def get_unique_lines():
    """only return the top hit for each query"""
    curr_dir=os.getcwd()
    for infile in glob.glob(os.path.join(curr_dir, "*.filtered")):
        names = get_seq_name(infile)
        outfile = open("%s.filtered.unique" % names, "w")
        d = {}
        input = file(infile)
        outdata = [ ]
        for line in input:
            unique = line.split("\t",1)[0]
            if unique not in d:
                d[unique] = 1
                outfile.write(line)
                outdata.append(line)
        outfile.close()
    return outdata

def blast_against_self_blastn(blast_type, genes_pep, genes_nt, output, filter, processors):
    devnull = open('/dev/null', 'w')
    if "F" in filter:
        my_seg = "yes"
    else:
        my_seg = "no"
    cmd = ["%s" % blast_type,
           "-task", blast_type,
           "-query", genes_pep,
           "-db", genes_nt,
           "-num_threads", str(processors),
           "-evalue", "0.1",
           "-outfmt", "6",
           "-dust", str(my_seg),
           "-out", output]
    subprocess.call(cmd, stdout=devnull, stderr=devnull)
    #subprocess.run(cmd,stdout=devnull,stderr=devnull)

def blast_against_self_tblastn(blast_type, genes_nt, genes_pep, output,processors, filter):
    devnull = open('/dev/null', 'w')
    if "F" in filter:
        my_seg = "no"
    else:
        my_seg = "yes"
    cmd = ["%s" % blast_type,
           "-query", genes_pep,
           "-db", genes_nt,
           "-num_threads", str(processors),
           "-seg", my_seg,
           "-comp_based_stats", "F",
           "-evalue", "0.1",
           "-outfmt", "6",
           "-out", output]
    try:
        subprocess.call(cmd, stdout=devnull, stderr=devnull)
    except:
        print("error!")

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

def translate_genes(genes,outfile,min_len):
    """translate nucleotide into peptide with BioPython"""
    infile = open(genes, "rU")
    output = []
    output_handle = open(outfile, "w")
    too_short = []
    for record in SeqIO.parse(infile, "fasta"):
        try:
            min_pep_len=int(min_len)
            """Should I trim these sequences back to be multiples of 3?"""
            if (len(record.seq)/3.0).is_integer():
                pep_seq=record.seq.translate(to_stop=True, table=11)
            elif ((len(record.seq)-1)/3.0).is_integer():
                pep_seq=record.seq[:-1].translate(to_stop=True, table=11)
            elif ((len(record.seq)-2)/3.0).is_integer():
                pep_seq=record.seq[:-2].translate(to_stop=True, table=11)
            elif ((len(record.seq)-3)/3.0).is_integer():
                pep_seq=record.seq[:-3].translate(to_stop=True, table=11)
            else:
                print("Sequence of odd length found and couldn't be trimmed")
            if len(pep_seq)>=min_pep_len:
                output_handle.write(">"+record.id+"\n")
                output_handle.write("".join(pep_seq)+"\n")
                output.append(pep_seq)
            else:
                too_short.append(record.id)
        except:
            raise TypeError("odd characters observed in sequence %s" % record.id)
    infile.close()
    output_handle.close()
    for record in output:
        return str(record)
    if len(too_short)>0:
        logging.logPrint("The following sequences were too short and will not be processed: %s" % "\n".join(too_short))

rec=1

def autoIncrement():
    global rec
    pStart = 1
    pInterval = 1
    if rec == 0:
        rec = pStart
    else:
        rec += pInterval
        return rec

def prune_matrix(matrix, group1, group2):
    """prune out genomes of interest from a BSR matrix.
    Not done efficiently, but appears to work"""
    in_matrix = open(matrix, "U")
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
    group1_out.write("\t"+"\t"+"\t".join(fields))
    for line in in_matrix:
        fields = line.split()
        name = fields[0]
        deque((list.pop(fields, i) for i in sorted(group1_idx, reverse=True)), maxlen=0)
    group1_out.write("".join(name)+"\t"+"\t".join(fields))
    in_matrix = open(matrix, "U")
    firstLine = in_matrix.readline()
    fields = firstLine.split()
    fields.insert(0, "cluster")
    group2_idx = [ ]
    for x in fields:
        if x not in group2_ids: group2_idx.append(fields.index(x))
    deque((list.pop(fields, i) for i in sorted(group2_idx, reverse=True)), maxlen=0)
    group2_out.write("\t"+"\t".join(fields))
    for line in in_matrix:
        fields = line.split()
        name = fields[0]
        deque((list.pop(fields, i) for i in sorted(group2_idx, reverse=True)), maxlen=0)
        group2_out.write("".join(name)+"\t"+"\t".join(fields))
    return group1_ids, group2_ids, group1_idx, group2_idx
    in_matrix.close()

def compare_values(pruned_1,pruned_2,upper,lower):
    group1 = open(pruned_1, "U")
    group2 = open(pruned_2, "U")
    group1_out = open("group1_out.txt", "w")
    group2_out = open("group2_out.txt", "w")
    group1_presents=[ ]
    group2_presents=[ ]
    group1_mean = [ ]
    next(group1)
    for line in group1:
        fields = line.split()
        presents = [ ]
        homolog = [ ]
        ints=map(float, fields[1:])
        mean = float(sum(ints)/len(ints))
        group1_mean.append(mean)
        for x in ints:
            if float(x)>=float(upper): presents.append(x)
            if float(x)>=float(upper): group1_presents.append(x)
            if float(x)>=float(lower): homolog.append(x)
        group1_out.write(str(fields[0])+"\t"+str(mean)+"\t"+str(len(presents))+"\t"+str(len(fields[1:]))+"\t"+str(len(homolog)))
    next(group2)
    for line in group2:
        fields = line.split()
        presents = [ ]
        homolog = [ ]
        ints=map(float, fields[1:])
        mean = float(sum(ints)/len(ints))
        for x in ints:
            if float(x)>=float(upper): presents.append(x)
            if float(x)>=float(upper): group2_presents.append(x)
            if float(x)>=float(lower): homolog.append(x)
        group2_out.write(str(mean)+"\t"+str(len(presents))+"\t"+str(len(fields[1:]))+"\t"+str(len(homolog)))
    return group1_presents, group2_presents, group1_mean
    group1.close()
    group2.close()
    group1_out.close()
    group2_out.close()

def find_uniques(combined,fasta):
    infile = open(combined, "U")
    group1_unique_ids = [ ]
    seqrecords=[ ]
    testids = [ ]
    for line in infile:
        fields=line.split()
        if int(fields[2])/int(fields[3])==1 and int(fields[8])==0:
            group1_unique_ids.append(fields[0])
    for record in SeqIO.parse(fasta, "fasta"):
        if record.id in group1_unique_ids:
            seqrecords.append(record)
            testids.append(record.id)
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
    return group1_unique_ids, group2_unique_ids, testids

def filter_genomes(genomes, in_matrix):
    in_matrix = open(in_matrix, "rU")
    firstLine = in_matrix.readline()
    first_fields = firstLine.split()
    all_genomes=first_fields
    genomes_file = open(genomes, "r").read().splitlines()
    genomes_file = [x.strip(' ') for x in genomes_file]
    to_keep = [ ]
    for x in all_genomes:
        if x in genomes_file:
            to_keep.append(all_genomes.index(x))
    return to_keep
    in_matrix.close()

def filter_matrix(to_keep, in_matrix, prefix):
    matrix = open(in_matrix, "rU")
    outfile = open("%s_genomes.matrix" % prefix, "w")
    outdata = [ ]
    to_remove = [x+1 for x in to_keep]
    firstLine = matrix.readline()
    first_fields = firstLine.split()
    deque((list.pop(first_fields, i) for i in sorted(to_keep, reverse=True)), maxlen=0)
    outdata.append(first_fields)
    first_fields.insert(0,"")
    outfile.write("\t".join(first_fields))
    for line in matrix:
        fields = line.split()
        deque((list.pop(fields, i) for i in sorted(to_remove, reverse=True)), maxlen=0)
        outfile.write("\t".join(fields))
        outdata.append(fields)
    outfile.close()
    return outdata

def get_core_gene_stats(matrix, threshold, lower):
    in_matrix=open(matrix, "U")
    outfile = open("core_gene_ids.txt", "w")
    singletons = open("unique_gene_ids.txt", "w")
    firstLine = in_matrix.readline()
    positives = [ ]
    singles = [ ]
    for line in in_matrix:
        fields = line.split()
        totals = len(fields[1:])
        presents = [ ]
        uniques = [ ]
        try:
            for x in fields[1:]:
                if float(x)>=float(threshold):
                    presents.append(fields[0])
                if float(x)>=float(lower):
                    uniques.append(fields[0])
            if int(len(presents))/int(totals)>=1:
                positives.append(fields[0])
            if int(len(uniques))==1:
                singles.append(fields[0])
        except:
            raise TypeError("problem in input file found")

    print("# of conserved genes = %s" % len(positives))
    print("# of unique genes = %s" % len(singles))
    ratio = int(len(singles))/int(totals)
    outfile.write("\n".join(positives)+"\n")
    singletons.write("\n".join(singles)+"\n")
    print("# of unique genes per genome = %s" % ratio)
    in_matrix.close()
    outfile.close()
    singletons.close()
    return len(positives), len(singles)

def get_frequencies(matrix, threshold):
    in_matrix=open(matrix, "U")
    firstLine = in_matrix.readline()
    outfile = open("frequency_data.txt", "w")
    my_dict = {}
    out_data = [ ]
    all = [ ]
    for line in in_matrix:
        presents = [ ]
        tempo = [ ]
        fields = line.split()
        try:
            for x in fields[1:]:
                if float(x)>=float(threshold):
                    presents.append(fields[0])
        except:
            raise TypeError("problem found with input file")
        tempo.append(len(presents))
        tempo.append("1")
        all.append(tempo)
    for x, y in all:
        try:
            my_dict[x].append(y)
        except KeyError:
            my_dict[x]=[y]
    outfile.write("Frequency distribution:\n")
    for k,v in my_dict.iteritems():
        outfile.write(str(k)+"\t"+str(len(v))+"\n")
        out_data.append(k)
        out_data.append(len(v))
    in_matrix.close()
    outfile.close()
    return out_data

def find_dups_dev(ref_scores, length, max_plog, min_hlog, clusters, processors):
    curr_dir=os.getcwd()
    my_dict_o = {}
    dup_dict = {}
    paralogs = [ ]
    duplicate_file = open("duplicate_ids.txt", "w")
    paralog_file = open("paralog_ids.txt", "w")
    ref_file = open("dup_refs.txt", "w")
    genome_specific_list_of_lists = []
    target_list = []
    files = os.listdir(curr_dir)
    files_and_temp_names = [(str(idx), os.path.join(curr_dir, f))
                            for idx, f in enumerate(files)]
    def _perform_workflow(data):
        tn, f = data
        if "_blast.out" in f:
            genome_specific_dict = {}
            name = get_seq_name(f)
            reduced_name = name.replace(".fasta.new_blast.out","")
            genome_specific_dict.update({"ID":reduced_name})
            outfile = open("%s.counts.txt" % reduced_name, "w")
            try:
                for line in open(f, "U"):
                    fields = line.split()
                    if fields[0] not in ref_scores: pass
                    elif float(fields[2])>=int(min_hlog) and (float(fields[11])/float(ref_scores.get(fields[0])))>=float(length):
                        try:
                            my_dict_o[fields[0]].append(fields[11])
                            genome_specific_dict[fields[0]].append(fields[11])
                        except KeyError:
                            my_dict_o[fields[0]] = [fields[11]]
                            genome_specific_dict[fields[0]] = [fields[11]]
                    else:
                        continue
            except:
                raise TypeError("problem parsing %s" % f)
            new_dict = {}
            for k,v in genome_specific_dict.iteritems():
                for cluster in clusters:
                    if k == "ID":
                        pass
                    elif k == cluster:
                        try:
                            new_dict.update({k:len(v)})
                        except:
                            new_dict.update({k:"0"})
            for cluster in clusters:
                if cluster not in genome_specific_dict:
                    new_dict.update({cluster:"0"})
            od = collections.OrderedDict(sorted(new_dict.items()))
            ids = OrderedDict({"ID":reduced_name})
            both =OrderedDict(list(ids.items())+list(new_dict.items()))
            for k,v in both.iteritems():
                outfile.write(str(v)+"\n")
                if k in target_list:
                    pass
                else:
                    target_list.append(k)
            outfile.close()
    results = set(p_func.pmap(_perform_workflow,
                              files_and_temp_names,
                              num_workers=processors))
    ref_file.write("\n".join(target_list)+"\n")
    ref_file.close()
    """known issue - if gene id is Capital and before "I", there can be a shuffling of IDs
    I need to sort the dictionary and keep the first item constant as ID"""
    try:
        generate_dup_matrix()
        os.system("paste dup_refs.txt dup_values > dup_matrix.txt")
    except:
        print("problem generating duplicate matrix, but we'll continue")
    for k,v in my_dict_o.iteritems():
        if int(len(v))>=2:
            dup_dict.update({k:v})
    for k,v in dup_dict.iteritems():
        max_value = max(v)
        for x in v:
            if float(x)/float(max_value)<=max_plog:
                paralogs.append(k)
            else:
                continue
    for k, v in dup_dict.iteritems():
        duplicate_file.write(k+"\n")
    nr=[x for i, x in enumerate(paralogs) if x not in paralogs[i+1:]]
    paralog_file.write("\n".join(nr)+"\n")
    duplicate_file.close()
    paralog_file.close()
    return nr, dup_dict

def filter_paralogs(matrix, ids):
    in_matrix = open(matrix, "U")
    outfile = open("bsr_matrix_values_filtered.txt", "w")
    outdata = [ ]
    genomes_file = open(ids, "rU").read().splitlines()
    firstLine = in_matrix.readline()
    outfile.write(firstLine)
    for line in in_matrix:
        fields = line.split()
        if fields[0] not in genomes_file:
            outfile.write(line)
            outdata.append(fields[0])
        else:
            pass
    return outdata
    in_matrix.close()
    outfile.close()

def filter_variome(matrix, threshold, step):
    in_matrix = open(matrix, "U")
    outfile = open("variome_BSR_matrix", "w")
    firstLine = in_matrix.readline()
    outdata = [ ]
    outfile.write(firstLine)
    for line in in_matrix:
        fields = line.split()
        totals = len(fields[1:])
        presents = [ ]
        for x in fields[1:]:
            try:
                if float(x)>=float(threshold):
                    presents.append(fields[0])
            except:
                raise TypeError("problem in input file observed")
        if int(len(presents))<(totals-int(step)):
            outdata.append(fields[0])
            outfile.write(line)
    in_matrix.close()
    outfile.close()
    return outdata

def run_usearch(id):
    rec=1
    curr_dir=os.getcwd()
    devnull = open("/dev/null", "w")
    for infile in glob.glob(os.path.join(curr_dir, "x*")):
        cmd = ["usearch",
           "-cluster_fast", "%s" % infile,
           "-id", str(id),
           "-uc", "results.uc",
           "-centroids", "%s.usearch.out" % str(autoIncrement())]
        subprocess.call(cmd,stdout=devnull,stderr=devnull)
    devnull.close()

def filter_scaffolds(in_fasta):
    """If an N is present in any scaffold, the entire contig will
    be entire filtered, probably too harsh"""
    infile = open(in_fasta, "U")
    outrecords = [ ]
    for record in SeqIO.parse(infile, "fasta"):
        if "N" not in record.seq:
            outrecords.append(record)
    output_handle = open("tmp.out", "w")
    if int(len(outrecords))==0:
        print("no usable fasta records were found or all contain scaffolds")
        sys.exit()
    SeqIO.write(outrecords, output_handle, "fasta")
    output_handle.close()

def uclust_sort(usearch):
    """sort with Usearch. Updated to V6"""
    devnull = open("/dev/null", "w")
    cmd = ["%s" % usearch,
           "-sortbylength", "all_gene_seqs.out",
           "-output", "tmp_sorted.txt"]
    subprocess.call(cmd,stdout=devnull,stderr=devnull)
    devnull.close()

def process_pangenome(matrix, upper, lower, iterations, type, prefix):
    my_matrix = open(matrix, "U")
    if "acc" in type:
        acc_outfile = open("%s_accumulation_replicates.txt" % prefix, "w")
    elif "uni" in type:
        uni_outfile = open("%s_uniques_replicates.txt" % prefix, "w")
    elif "core" in type:
        core_outfile = open("%s_core_replicates.txt" % prefix, "w")
    else:
        acc_outfile = open("%s_accumulation_replicates.txt" % prefix, "w")
        uni_outfile = open("%s_uniques_replicates.txt" % prefix, "w")
        core_outfile = open("%s_core_replicates.txt" % prefix, "w")
    firstLine = my_matrix.readline()
    first_fields = firstLine.split()
    genomes = len(first_fields)
    indexes = []
    for x in first_fields:
        indexes.append(first_fields.index(x)+1)
    my_matrix.close()
    acc_dict = {}
    core_dict = {}
    uni_dict = {}
    for j in range(1,iterations+1):
        for i in range(1,genomes+1):
            positives_acc = []
            positives_core = []
            positives_unis = []
            outseqs=random.sample(set(indexes), int(i))
            with open(matrix, "U") as f:
                next(f)
                for line in f:
                    fields = line.split()
                    positive_lines_acc=[]
                    positive_lines_core=[]
                    positive_lines_unis=[]
                    for outseq in outseqs:
                        if type == "acc" or type == "all":
                            if float(fields[outseq])>=float(upper):
                                positive_lines_acc.append("1")
                        if type == "core" or type == "all":
                            if float(fields[outseq])>=float(upper):
                                positive_lines_core.append("1")
                        if type == "uni" or type == "all":
                            """this was changed from lower to upper"""
                            if float(fields[outseq])>=float(lower) and float(fields[outseq])>=float(upper):
                                positive_lines_unis.append("1")
                    if len(positive_lines_acc)>=1:
                        positives_acc.append("1")
                    if len(positive_lines_core)==len(outseqs):
                        positives_core.append("1")
                    if int(len(positive_lines_unis))==1:
                        positives_unis.append("1")
            try:
                acc_dict[i].append(len(positives_acc))
            except KeyError:
                acc_dict[i] = [len(positives_acc)]
            try:
                core_dict[i].append(len(positives_core))
            except KeyError:
                core_dict[i] = [len(positives_core)]
            try:
                uni_dict[i].append(len(positives_unis))
            except KeyError:
                uni_dict[i] = [len(positives_unis)]
    try:
        sorted_acc_dict = collections.OrderedDict(sorted(acc_dict.items()))
        sorted_uni_dict = collections.OrderedDict(sorted(uni_dict.items()))
        sorted_core_dict = collections.OrderedDict(sorted(core_dict.items()))
    except:
        pass
    test_accums = []
    test_uniques = []
    test_cores = []
    if type == "acc" or type == "all":
        print("accumulation means")
        for k,v in sorted_acc_dict.iteritems():
            test_accums.append(v)
            print(k, sum(v)/len(v))
            for z in v:
                acc_outfile.write(str(k)+"\t"+str(z)+"\n")
    if type == "uni" or type == "all":
        print("unique means")
        for k,v in sorted_uni_dict.iteritems():
            test_uniques.append(v)
            print(k, (sum(v)/len(v))/int(k))
            for z in v:
                uni_outfile.write(str(k)+"\t"+str(int(z)/int(k))+"\n")
    if type == "core" or type == "all":
        print("core means")
        for k,v in sorted_core_dict.iteritems():
            test_cores.append(v)
            print(k, sum(v)/len(v))
            for z in v:
                core_outfile.write(str(k)+"\t"+str(z)+"\n")
    try:
        acc_outfile.close()
        uni_outfile.close()
        core_outfile.close()
    except:
        pass
    return test_accums, test_uniques, test_cores

def bsr_to_pangp(matrix, lower):
    my_matrix = open(matrix, "U")
    outfile = open("panGP_matrix.txt","w")
    firstLine = my_matrix.readline()
    outfile.write(firstLine)
    for line in my_matrix:
        new_fields = [ ]
        fields = line.split()
        new_fields.append(fields[0])
        for x in fields[1:]:
            if float(x)>=float(lower):
                new_fields.append("1")
            else:
                new_fields.append("-")
        outfile.write("\t".join(new_fields))
    my_matrix.close()
    outfile.close()
    return new_fields

def transpose_matrix(matrix):
    out_matrix = open("tmp.matrix", "w")
    reduced = [ ]
    for line in open(matrix, "U"):
        newline=line.strip("\n")
        fields = newline.split("\t")
        reduced.append(fields)
    test=map(list, zip(*reduced))
    for x in test:
        out_matrix.write("\t".join(x))
    out_matrix.close()

def reorder_matrix(in_matrix, names):
    my_matrix = open(in_matrix, "U")
    outfile = open("reordered_matrix.txt", "w")
    firstLine = my_matrix.readline()
    outfile.write(firstLine)
    my_matrix.close()
    for name in names:
         for line in open(in_matrix, "U"):
            newline = line.strip("\n")
            fields = newline.split("\t")
            if name == fields[0]:
                outfile.write(line)
    my_matrix.close()
    outfile.close()

def parse_tree(tree):
    names = []
    mytree = Phylo.read(tree, 'newick')
    for clade in mytree.find_clades():
        if clade.name:
            names.append(clade.name)
    return names

def blat_against_self(query,reference,output,processors):
    subprocess.check_call("blat -out=blast8 -minIdentity=75 %s %s %s > /dev/null 2>&1" % (reference,query,output), shell=True)

def blat_against_each_genome(dir_path,database,processors):
    """BLAT all genes against each genome"""
    curr_dir=os.getcwd()
    files = os.listdir(curr_dir)
    files_and_temp_names = [(str(idx), os.path.join(curr_dir, f))
                            for idx, f in enumerate(files)]
    def _perform_workflow(data):
        tn, f = data
        if ".fasta.new" in f:
            try:
                subprocess.check_call("blat -out=blast8 -minIdentity=75 %s %s %s_blast.out > /dev/null 2>&1" % (f,database,f), shell=True)
            except:
                print("genomes %s cannot be used" % f)

    results = set(p_func.pmap(_perform_workflow,
                              files_and_temp_names,
                              num_workers=processors))

def make_table_dev(infile, test, clusters):
    """make the BSR matrix table"""
    values = [ ]
    names = [ ]
    outdata = [ ]
    name=[ ]
    out=get_seq_name(infile)
    name.append(out)
    reduced=[ ]
    """remove the junk at the end of the file"""
    for x in name:reduced.append(x.replace('.fasta.new_blast.out.filtered.filtered.unique',''))
    names.append(reduced)
    my_dict={}
    my_file=open(infile, "rU")
    """make a dictionary of all clusters and values"""
    try:
        for line in my_file:
            fields=line.split()
            my_dict.update({fields[0]:fields[1]})
    except:
        raise TypeError("abnormal number of fields")
    my_file.close()
    """add in values, including any potentially missing ones"""
    for x in clusters:
        if x not in my_dict.keys():my_dict.update({x:0})
    for x in reduced:
        values.append(x)
    """sort keys to get the same order between samples"""
    od = collections.OrderedDict(sorted(my_dict.items()))
    values_2 = od.values()
    values_3 = values+values_2
    if "T" in test:
        myout=[x for i, x in enumerate(outdata) if x not in outdata[i+1:]]
        return sorted(outdata)
    else:
        pass
    return names, values_3

def create_bsr_matrix_dev(master_list):
    new_matrix = open("bsr_matrix", "w")
    test = map(list, zip(*master_list))
    for x in test:
        y = map(str, x)
        new_matrix.write("\t".join(y)+"\n")
    new_matrix.close()

def new_loop(to_iterate, processors, clusters, debug):
    names = []
    table_list = []
    def _perform_workflow(data):
        tn, f = data
        name,values=make_table_dev(f, "F", clusters)
        names.append(name)
        table_list.append(values)
        if debug == "T":
            logging.logPrint("sample %s processed" % f)
        else:
            pass
    set(p_func.pmap(_perform_workflow,
                    to_iterate,
                    num_workers=processors))
    return names,table_list

def run_vsearch(id, processors):
    devnull = open("/dev/null", "w")
    cmd = ["vsearch",
           "-cluster_fast", "all_gene_seqs.out",
           "-id", str(id),
           "-uc", "results.uc",
           "-threads", "%s" % processors,
           "-centroids", "vsearch.out"]

    subprocess.call(cmd,stdout=devnull,stderr=devnull)
    devnull.close()

def process_genbank_files(directory):
    genbank_hits = []
    for infile in glob.glob(os.path.join(directory, "*.gbk")):
        name = get_seq_name(infile)
        reduced = name.replace(".gbk","")
        genbank_hits.append(name)
        record = SeqIO.read(infile, "genbank")
        output_handle = open("%s.locus_tags.fasta" % reduced, "w")
        count = 0
        for feature in record.features:
            if feature.type == "gene":
                count = count + 1
                try:
                    feature_name = feature.qualifiers["locus_tag"]
                    feature_seq = feature.extract(record.seq)
                    output_handle.write(">" + "".join(feature_name) + "\n" + str(feature_seq) + "\n")
                except:
                    print("problem extracting locus tag: %s" % "".join(feature_name))
        output_handle.close()
    return genbank_hits

def test_duplicate_header_ids(fasta_file):
    IDs = []
    for line in open(fasta_file):
        if line.startswith(">"):
            fields = line.split()
            clean = fields[0].replace(">","")
            IDs.append(clean)
        else:
            pass
    nr=[x for i, x in enumerate(IDs) if x not in IDs[i+1:]]
    if len(IDs) == len(nr):
        return "True"
    else:
        return "False"

def split_files(fasta_file):
    """This next section removes line wraps, so I can
    split the file without interrupting a gene"""
    from Bio.SeqIO.FastaIO import FastaWriter
    output_handle = open("nowrap.fasta", "w")
    seqrecords=[ ]
    writer = FastaWriter(output_handle, wrap=0)
    for record in SeqIO.parse(open(fasta_file), "fasta"):
        seqrecords.append(record)
    writer.write_file(seqrecords)
    output_handle.close()
    """I can always make the number of lines an alterable field"""
    subprocess.check_call("split -l 200000 nowrap.fasta", shell=True)

def generate_dup_matrix():
    curr_dir=os.getcwd()
    dup_names = []
    outfile = open("dup_values", "w")
    for infile in glob.glob(os.path.join(curr_dir, '*.counts.txt')):
        genome_fields = []
        for line in open(infile, "rU"):
            newline = line.strip()
            fields = newline.split()
            for field in fields:
                genome_fields.append(field)
        dup_names.append(genome_fields)
    test=map(list, zip(*dup_names))
    for alist in test:
        outfile.write("\t".join(alist)+"\n")
    outfile.close()
