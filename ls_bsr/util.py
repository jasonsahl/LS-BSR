#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
import sys
import time
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
#import igs_logging as logging
import errno
import threading
import types
from collections import deque,OrderedDict
import collections

def mp_shell(func, params, numProc):
    from multiprocessing import Pool
    p = Pool(numProc)
    out = p.map(func, params)
    p.terminate()
    return out

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

def divide_values(file, ref_scores):
    """divide each BSR value in a row by that row's maximum value"""
    errors = []
    outdata = []
    with open(file) as infile:
        firstLine = infile.readline()
        FL_F=firstLine.split()
        outfile = open("BSR_matrix_values.txt", "w")
        outfile.write('\t'.join([str(item) for item in FL_F])+"\n")
        for line in infile:
            fields=line.split()
            all_fields=list(fields)
            try:
                fields=list(map(float, all_fields[1:]))
            except:
                raise TypeError("abnormal number of fields observed")
            values= []
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
        outfile.close()
    if len(errors)>0:
        nr=[x for i, x in enumerate(errors) if x not in errors[i+1:]]
        logPrint("The following genes had no hits in datasets or are too short, values changed to 0, check names and output:%s" % "\n".join(nr))
    return outdata

def rename_fasta_header(fasta_in, fasta_out):
    """this is used for renaming the output,
    in the off chance that there are duplicate
    names for separate peptides"""
    rec=1
    handle = open(fasta_out, "w")
    outdata = []
    with open(fasta_in) as infile:
        for record in SeqIO.parse(infile, "fasta"):
            try:
                outdata.append(">"+"centroid"+"_"+record.id)
                handle.write(">centroid_"+str(autoIncrement())+"\n")
                handle.write(str(record.seq)+"\n")
            except:
                raise TypeError("problem with input sequence encountered")
    handle.close()
    return outdata

def uclust_cluster(id, data_type):
    devnull = open("/dev/null", "w")
    if "nt" in data_type:
        output_name = "consensus.fasta"
    else:
        output_name = "consensus.pep"
    """cluster with Uclust.  Updated to V6"""
    cmd = ["usearch",
           "-cluster_fast", "all_sorted.txt",
           "-id", str(id),
           "-uc", "results.uc",
           "-centroids", str(output_name)]
    subprocess.call(cmd, stderr=devnull, stdout=devnull)

def get_seq_name(in_fasta):
    """used for renaming the sequences"""
    return os.path.basename(in_fasta)

def filter_seqs(input_pep,output_pep):
    """filter out short sequences from a multifasta.
    Will hopefully speed up the process without losing
    important information"""
    long_sequences = []
    outfile = open(output_pep, "w")
    outdata = []
    with open(input_pep) as infile:
        for record in SeqIO.parse(infile, "fasta"):
            if len(record.seq) >= int(50):
                long_sequences.append(record)
                outdata.append(len(record.seq))
    SeqIO.write(long_sequences, outfile, "fasta")
    outfile.close()
    return outdata

def parse_blast_report_dev(test,processors):
    """parse out only the unqiue names and bit score from the blast report"""
    curr_dir=os.getcwd()
    files_and_temp_names = []
    for infile in glob.glob(os.path.join(curr_dir, "*_blast.out")):
        files_and_temp_names.append([infile, test])
    outdata = mp_shell(_perform_workflow_pbr, files_and_temp_names, processors)
    if "true" in test:
        # mp_shell will return a list of lists. This will flatten it into a single list
        return outdata
    return

def _perform_workflow_pbr(data):
    infile = data[0]
    test = data[1]
    outdata = []
    order = []
    names = get_seq_name(infile)
    outfile = open("%s.filtered.unique" % names, "w")
    uniques = {}
    with open(data[0]) as infile:
        for line in infile:
            try:
                fields = line.split()
                # Keep track of the largest value of fields[0]
                if fields[0] not in uniques:
                    uniques[fields[0]] = fields[11].strip("\n")
                    order.append(fields[0])
                else:
                    if float(fields[11]) > float(uniques[fields[0]]):
                        uniques[fields[0]] = fields[11].strip("\n")
            except IndexError:
                raise TypeError("Malformed blast line found in %s" % infile)
    for item in order:
        if "true" in test:
            outdata.append(item)
            outdata.append(uniques[item])
        outfile.write(item + "\t" + uniques[item] + "\n")
    outfile.close()
    if "true" in test:
        return outdata

def blast_against_self_blastn(blast_type, algorithm, genes_pep, genes_nt, output, filter, processors):
    devnull = open('/dev/null', 'w')
    if "F" in filter:
        my_seg = "yes"
    else:
        my_seg = "no"
    cmd = ["%s" % blast_type,
           "-task", algorithm,
           "-query", genes_pep,
           "-db", genes_nt,
           "-num_threads", str(processors),
           "-evalue", "0.1",
           "-outfmt", "6",
           "-dust", str(my_seg),
           "-out", output]
    subprocess.call(cmd, stdout=devnull, stderr=devnull)

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

def parse_self_blast(infile):
    my_dict={}
    with open(infile) as my_blast:
        for line in my_blast:
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
    output = []
    output_handle = open(outfile, "w")
    too_short = []
    with open(genes) as infile:
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
    output_handle.close()
    for record in output:
        return str(record)
    if len(too_short)>0:
        logPrint("The following sequences were too short and will not be processed: %s" % "\n".join(too_short))

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
    group1_ids = []
    group2_ids = []
    group1_out = open("group1_pruned.txt", "w")
    group2_out = open("group2_pruned.txt", "w")
    with open(group1) as file_1:
        for line in file_1:
            line.strip()
            group1_ids.append(line)
    with open(group2) as file_2:
        for line in file_2:
            line.strip()
            group2_ids.append(line)
    with open(matrix) as in_matrix:
        """This is just to get the headers correct"""
        firstLine = in_matrix.readline()
        fields = firstLine.split()
        fields.insert(0, "cluster")
        group1_ids = map(lambda s: s.strip(), group1_ids)
        group2_ids = map(lambda s: s.strip(), group2_ids)
        group1_idx = []
        group1_ids_list = list(group1_ids)
        group2_ids_list = list(group2_ids)
        for x in fields:
            """These are the genomes that I want to prune out"""
            if x not in group1_ids_list: group1_idx.append(fields.index(x))
        deque((list.pop(fields, i) for i in sorted(group1_idx, reverse=True)), maxlen=0)
        group1_out.write("\t"+"\t"+"\t".join(fields)+"\n")
        for line in in_matrix:
            fields = line.split()
            name = fields[0]
            if name == "cluster":
                pass
            else:
                try:
                    deque((list.pop(fields, i) for i in sorted(group1_idx, reverse=True)), maxlen=0)
                    group1_out.write("".join(name)+"\t"+"\t".join(fields)+"\n")
                except:
                    pass
    group1_out.close()
    with open(matrix) as in_matrix:
        firstLine = in_matrix.readline()
        fields = firstLine.split()
        fields.insert(0, "cluster")
        group2_idx = []
        for x in fields:
            if x not in group2_ids_list: group2_idx.append(fields.index(x))
        deque((list.pop(fields, i) for i in sorted(group2_idx, reverse=True)), maxlen=0)
        group2_out.write("\t"+"\t".join(fields)+"\n")
        for line in in_matrix:
            fields = line.split()
            name = fields[0]
            if name == "cluster":
                pass
            else:
                try:
                    deque((list.pop(fields, i) for i in sorted(group2_idx, reverse=True)), maxlen=0)
                    group2_out.write("".join(name)+"\t"+"\t".join(fields)+"\n")
                except:
                    pass
    group2_out.close()
    return group1_ids_list, group2_ids_list, group1_idx, group2_idx


def compare_values(pruned_1,pruned_2,upper,lower):
    group1_out = open("group1_out.txt", "w")
    group2_out = open("group2_out.txt", "w")
    group1_presents=[]
    group2_presents=[]
    group1_mean = []
    with open(pruned_1) as group1:
        next(group1)
        for line in group1:
            fields = line.split()
            presents = []
            homolog = []
            if len(fields) == 0:
                pass
            else:
                ints=list(map(float, fields[1:]))
                num_ints = len(ints)
                sum_ints = sum(ints)
                mean = float(sum_ints/num_ints)
                group1_mean.append(mean)
                for x in ints:
                    if float(x)>=float(upper): presents.append(x)
                    if float(x)>=float(upper): group1_presents.append(x)
                    if float(x)>=float(lower): homolog.append(x)
                group1_out.write(str(fields[0])+"\t"+str(mean)+"\t"+str(len(presents))+"\t"+str(len(fields[1:]))+"\t"+str(len(homolog))+"\n")
    with open(pruned_2) as group2:
        next(group2)
        for line in group2:
            fields = line.split()
            presents = []
            homolog = []
            if len(fields) == 0:
                pass
            else:
                ints=list(map(float, fields[1:]))
                num_ints = len(ints)
                sum_ints = sum(ints)
                mean = float(sum_ints/num_ints)
                for x in ints:
                    if float(x)>=float(upper): presents.append(x)
                    if float(x)>=float(upper): group2_presents.append(x)
                    if float(x)>=float(lower): homolog.append(x)
                group2_out.write(str(mean)+"\t"+str(len(presents))+"\t"+str(len(fields[1:]))+"\t"+str(len(homolog))+"\n")
    group1_out.close()
    group2_out.close()
    return group1_presents, group2_presents, group1_mean

def find_uniques(combined,fasta):
    group1_unique_ids = []
    seqrecords=[]
    testids = []
    with open(combined) as infile:
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
    group2_unique_ids = []
    seqrecords2 = []
    with open(combined) as infile:
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
    to_keep = []
    with open(in_matrix) as matrix:
        firstLine = matrix.readline()
        first_fields = firstLine.split()
        all_genomes=first_fields
        genomes_list = []
        with open(genomes) as genomes_file:
            for line in genomes_file:
                newline = line.strip()
                genomes_list.append(newline)
        genomes_file = [x.strip(' ') for x in genomes_list]
        for x in all_genomes:
            if x in genomes_file:
                to_keep.append(all_genomes.index(x))
    return to_keep

def filter_matrix(to_keep, in_matrix, prefix):
    with open(in_matrix) as matrix:
        outfile = open("%s_genomes.matrix" % prefix, "w")
        outdata = []
        to_remove = [x+1 for x in to_keep]
        firstLine = matrix.readline()
        first_fields = firstLine.split()
        deque((list.pop(first_fields, i) for i in sorted(to_keep, reverse=True)), maxlen=0)
        outdata.append(first_fields)
        first_fields.insert(0,"")
        outfile.write("\t".join(first_fields)+"\n")
        for line in matrix:
            fields = line.split()
            deque((list.pop(fields, i) for i in sorted(to_remove, reverse=True)), maxlen=0)
            outfile.write("\t".join(fields)+"\n")
            outdata.append(fields)
    outfile.close()
    return outdata

def get_core_gene_stats(matrix, threshold, lower, missing):
    outfile = open("core_gene_ids.txt", "w")
    singletons = open("unique_gene_ids.txt", "w")
    positives = []
    singles = []
    with open(matrix) as in_matrix:
        firstLine = in_matrix.readline()
        for line in in_matrix:
            fields = line.split()
            totals = len(fields[1:])
            presents = []
            uniques = []
            try:
                for x in fields[1:]:
                    if float(x)>=float(threshold):
                        presents.append(fields[0])
                    if float(x)>=float(lower):
                        uniques.append(fields[0])
                if int(len(presents)+missing)/int(totals)>=1:
                    positives.append(fields[0])
                if int(len(uniques))==1:
                    singles.append(fields[0])
            except:
                raise TypeError("problem in input file found")
    print("# of conserved genes (>=0.8 BSR in all genomes) = %s" % len(positives))
    print("# of unique genes (>=0.8 BSR in only 1 genome, <0.4 in others) = %s" % len(singles))
    ratio = int(len(singles))/int(totals)
    outfile.write("\n".join(positives)+"\n")
    singletons.write("\n".join(singles)+"\n")
    print("# of unique genes per genome = %s" % ratio)
    outfile.close()
    singletons.close()
    return len(positives),len(singles)

def get_frequencies(matrix, threshold):
    import collections
    out_data = []
    all = []
    outfile = open("frequency_data.txt", "w")
    with open(matrix) as in_matrix:
        firstLine = in_matrix.readline()
        for line in in_matrix:
            """Presents contain those CDSs that are conserved"""
            presents = []
            tempo = []
            fields = line.split()
            try:
                for x in fields[1:]:
                    if float(x)>=float(threshold):
                        presents.append(fields[0])
            except:
                raise TypeError("problem found with input file")
            tempo.append(len(presents))
            all.append(tempo)
    """This now flattens the list of lists: There should never be zeros"""
    all_flat = [item for sublist in all for item in sublist]
    my_dict=collections.Counter(all_flat)
    outfile.write("Frequency distribution:\n")
    for k,v in my_dict.items():
        outfile.write(str(k)+"\t"+str(v)+"\n")
    outfile.close()
    return my_dict

def filter_paralogs(matrix, ids):
    outfile = open("bsr_matrix_values_filtered.txt", "w")
    my_list = []
    with open(ids) as genomes_file:
        for line in genomes_file:
            newline = line.strip()
            my_list.append(newline)
    num_to_filter = []
    with open(matrix) as in_matrix:
        firstLine = in_matrix.readline()
        outfile.write(firstLine)
        for line in in_matrix:
            fields = line.split()
            if fields[0] not in my_list:
                outfile.write(line)
            else:
                num_to_filter.append("1")
                pass
    num_filtered = len(num_to_filter)
    outfile.close()
    return num_filtered

def filter_variome(matrix, threshold, step):
    with open(matrix) as in_matrix:
        outfile = open("variome_BSR_matrix", "w")
        firstLine = in_matrix.readline()
        outdata = []
        outfile.write(firstLine)
        for line in in_matrix:
            fields = line.split()
            totals = len(fields[1:])
            presents = []
            for x in fields[1:]:
                try:
                    if float(x)>=float(threshold):
                        presents.append(fields[0])
                except:
                    raise TypeError("problem in input file observed")
            if int(len(presents))<(totals-int(step)):
                outdata.append(fields[0])
                outfile.write(line)
        outfile.close()
        return outdata

def filter_scaffolds_fun(in_fasta):
    """If an N is present in any scaffold, the entire contig will
    be entire filtered, probably too harsh"""
    with open(in_fasta) as infile:
        outrecords = []
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
    indexes = []
    with open(matrix) as my_matrix:
        firstLine = my_matrix.readline()
        first_fields = firstLine.split()
        genomes = len(first_fields)
        #indexes = []
        for x in first_fields:
            indexes.append(first_fields.index(x)+1)
    acc_dict = {}
    core_dict = {}
    uni_dict = {}
    #I could step through these by 2?
    for i in range(1,genomes+1):
        for j in range(1,iterations+1):
            positives_acc = []
            positives_core = []
            positives_unis = []
            """This selects the random set of genomes"""
            outseqs=random.sample(set(indexes), int(i))
            """Changing the tabs for all of these genomes"""
            with open(matrix) as f:
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
    #Changed spacing
    test_accums = []
    test_uniques = []
    test_cores = []
    if type == "acc" or type == "all":
        print("accumulation means")
        for k,v in sorted_acc_dict.items():
            test_accums.append(v)
            print(k, sum(v)/len(v))
            for z in v:
                acc_outfile.write(str(k)+"\t"+str(z)+"\n")
    if type == "uni" or type == "all":
        print("unique means")
        for k,v in sorted_uni_dict.items():
            test_uniques.append(v)
            print(k, (sum(v)/len(v))/int(k))
            for z in v:
                uni_outfile.write(str(k)+"\t"+str(int(z)/int(k))+"\n")
    if type == "core" or type == "all":
        print("core means")
        for k,v in sorted_core_dict.items():
            test_cores.append(v)
            print(k,sum(v)/len(v))
            for z in v:
                core_outfile.write(str(k)+"\t"+str(z)+"\n")
    if "acc" in type:
        acc_outfile.close()
    elif "uni" in type:
        uni_outfile.close()
    elif "core" in type:
        core_outfile.close()
    else:
        acc_outfile.close()
        uni_outfile.close()
        core_outfile.close()
    return test_accums, test_uniques, test_cores

def bsr_to_pangp(matrix, lower):
    with open(matrix) as my_matrix:
        outfile = open("panGP_matrix.txt","w")
        firstLine = my_matrix.readline()
        outfile.write(firstLine)
        for line in my_matrix:
            new_fields = []
            fields = line.split()
            new_fields.append(fields[0])
            for x in fields[1:]:
                if float(x)>=float(lower):
                    new_fields.append("1")
                else:
                    new_fields.append("-")
            outfile.write("\t".join(new_fields)+"\n")
        outfile.close()
        return new_fields

def transpose_matrix(matrix):
    out_matrix = open("tmp.matrix", "w")
    reduced = []
    with open(matrix) as infile:
        for line in infile:
            newline=line.strip("\n")
            fields = newline.split("\t")
            reduced.append(fields)
    test=map(list, zip(*reduced))
    for x in test:
        out_matrix.write("\t".join(x)+"\n")
    out_matrix.close()

def reorder_matrix(in_matrix, names):
    with open(in_matrix) as my_matrix:
        outfile = open("reordered_matrix.txt", "w")
        firstLine = my_matrix.readline()
        outfile.write(firstLine)
        my_matrix.close()
        for name in names:
            with open(in_matrix) as my_matrix:
                for line in my_matrix:
                    newline = line.strip("\n")
                    fields = newline.split()
                    if fields[0] in name:
                        outfile.write(line,)
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

def create_bsr_matrix_dev(master_list):
    new_matrix = open("bsr_matrix", "w")
    test = map(list,zip(*master_list))
    for x in test:
        y = map(str,x)
        new_matrix.write("\t".join(y)+"\n")
    new_matrix.close()

def run_vsearch(id, processors, infile):
    devnull = open("/dev/null", "w")
    cmd = ["vsearch",
           "-cluster_fast", infile,
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
        output_handle = open("%s.locus_tags.fasta" % reduced, "w")
        count = 0
        for record in SeqIO.parse(infile, "genbank"):
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
    with open(fasta_file) as infile:
        for line in infile:
            if line.startswith(">"):
                fields = line.split()
                clean = fields[0].replace(">","")
                IDs.append(clean)
            else:
                pass
    nr=[x for i, x in enumerate(IDs) if x not in IDs[i+1:]]
    my_dict = {}
    for id in IDs:
        try:
            my_dict[id].append("1")
        except KeyError:
            my_dict[id] = ["1"]
    dups = []
    for k,v in my_dict.items():
        if len(v)>1: dups.append(k)
    if len(dups)>0:
        print("Duplicate header IDs:")
        print("\n".join(dups))
    if len(IDs) == len(nr):
        return "True"
    else:
        return "False"

def split_files(fasta_file):
    """This next section removes line wraps, so I can
    split the file without interrupting a gene"""
    output_handle = open("nowrap.fasta", "w")
    with open(fasta_file) as infile:
        for record in SeqIO.parse(infile, "fasta"):
            output_handle.write(">"+str(record.id)+"\n")
            output_handle.write(str(record.seq)+"\n")
    output_handle.close()
    """I can always make the number of lines an alterable field"""
    subprocess.check_call("split -l 200000 nowrap.fasta", shell=True)

def generate_dup_matrix():
    curr_dir=os.getcwd()
    dup_names = []
    outfile = open("dup_values", "w")
    for infile in glob.glob(os.path.join(curr_dir, '*.counts.txt')):
        genome_fields = []
        with open(infile) as myfile:
            for line in myfile:
                newline = line.strip()
                fields = newline.split()
                for field in fields:
                    genome_fields.append(field)
        dup_names.append(genome_fields)
    test=map(list, zip(*dup_names))
    for alist in test:
        outfile.write("\t".join(alist)+"\n")
    outfile.close()

def _usearch_workflow(infile):
    devnull = open("/dev/null", "w")
    cmd = ["usearch",
           "-cluster_fast", "%s" % infile[0],
           "-id", str(infile[1]),
           "-sort", "length",
           "-uc", "results.uc",
           "-centroids", "%s.usearch.out" % infile[0]]
    subprocess.call(cmd,stdout=devnull,stderr=devnull)
    devnull.close()

def run_usearch_dev(id,processors):
    curr_dir=os.getcwd()
    # Put all files that start with 'x' in list
    files_and_temp_names = []
    for file in glob.glob(os.path.join(curr_dir, "x*")):
        files_and_temp_names.append([file,id])
    mp_shell(_usearch_workflow, files_and_temp_names, processors)

def _prodigal_workflow_def(data):
    tn, f = data
    subprocess.check_call("prodigal -i %s -d %s_genes.seqs -m -c -a %s_genes.pep > /dev/null 2>&1" % (f, f, f), shell=True)

def _prodigal_workflow_inter(data):
    tn,f = data
    name = f.replace(".fasta.new","")
    subprocess.check_call("prodigal -i %s -d %s_genes.seqs -a %s_genes.pep -f gff -m -c -o %s.prodigal > /dev/null 2>&1" % (f, f, f, name), shell=True)
    inverse_coding_regions("%s.prodigal" % name, name)
    parse_ranges_file(f,"%s.ranges" % name,name,test="false")

def predict_genes(fastadir, processors, intergenics):
    """simple gene prediction using Prodigal in order
    to find coding regions from a genome sequence"""
    os.chdir("%s" % fastadir)
    files = []
    for file in os.listdir(fastadir):
        if file.endswith(".fasta.new"):
            files.append(file)
    files_and_temp_names = [(str(idx), os.path.join(fastadir, f))
                            for idx, f in enumerate(files)]
    if intergenics == "F":
        mp_shell(_prodigal_workflow_def, files_and_temp_names, processors)
    else:
        mp_shell(_prodigal_workflow_inter, files_and_temp_names, processors)

def _perform_workflow_blat_genome(data):
    tn = data[0]
    f = data[1]
    database = data[2]
    if ".fasta.new" in f:
        try:
            subprocess.check_call("blat -out=blast8 -minIdentity=75 %s %s %s_blast.out > /dev/null 2>&1" % (f,database,f), shell=True)
        except:
            print("genomes %s cannot be used" % f)

def blat_against_each_genome_dev(database,processors):
    """BLAT all genes against each genome"""
    curr_dir=os.getcwd()
    files = []
    for file in os.listdir(curr_dir):
        if file.endswith(".fasta.new"):
            files.append(file)
    files_and_temp_names = []
    for idx,f in enumerate(files):
        files_and_temp_names.append([str(idx), os.path.join(curr_dir, f), database])
    mp_shell(_perform_workflow_blat_genome,files_and_temp_names,processors)

def _perform_workflow_tblastn(data):
    tn = data[0]
    f = data[1]
    my_seg = data[2]
    peptides = data[3]
    if f.endswith(".fasta.new"):
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
            devnull.close()
        except:
            print("genomes %s cannot be used" % f)

def blast_against_each_genome_tblastn_dev(processors, peptides, filter):
    """BLAST all peptides against each genome"""
    curr_dir=os.getcwd()
    files = []
    for file in os.listdir(curr_dir):
        if file.endswith(".fasta.new"):
            files.append(file)
    if "T" in filter:
        my_seg = "yes"
    else:
        my_seg = "no"
    files_and_temp_names = []
    for idx, f in enumerate(files):
        files_and_temp_names.append([str(idx), os.path.join(curr_dir, f), my_seg, peptides])
    mp_shell(_perform_workflow_tblastn, files_and_temp_names, processors)

def _perform_workflow_diamond(data):
    tn = data[0]
    f = data[1]
    peptides = data[2]
    name = f.replace(".new_genes.pep",".new")
    try:
        subprocess.check_call("diamond makedb --in %s -d %s > /dev/null 2>&1" % (f,name), shell=True)
    except:
        print("problem found in formatting annotation %s" % f)
    devnull = open('/dev/null', 'w')
    cmd = ["diamond",
           "blastp",
           "-p", "1",
           "-d", name,
           "-f", "6",
           "-q", peptides,
           "-o", "%s_blast.out" % name]
    subprocess.call(cmd, stdout=devnull, stderr=devnull)
    devnull.close()

def diamond_against_each_annotation(peptides,processors):
    curr_dir=os.getcwd()
    files_and_temp_names = []
    annotation_files = []
    for files in os.listdir(curr_dir):
        if "new_genes.pep" in files:
            annotation_files.append(files)
    for idx, f in enumerate(annotation_files):
        files_and_temp_names.append([str(idx), os.path.join(curr_dir, f), peptides])
    mp_shell(_perform_workflow_diamond, files_and_temp_names, processors)

def blastp_against_each_annotation(peptides,processors,filter):
    curr_dir=os.getcwd()
    files_and_temp_names = []
    annotation_files = []
    if "T" in filter:
        my_seg = "yes"
    else:
        my_seg = "no"
    for files in os.listdir(curr_dir):
        if "new_genes.pep" in files:
            annotation_files.append(files)
    for idx, f in enumerate(annotation_files):
        files_and_temp_names.append([str(idx), os.path.join(curr_dir, f), my_seg, peptides])
    mp_shell(_perform_workflow_blastp, files_and_temp_names, processors)

def _perform_workflow_blastp(data):
    tn = data[0]
    f = data[1]
    my_seg = data[2]
    peptides = data[3]
    """Makes the name consistent with other analyses"""
    name = f.replace(".new_genes.pep",".new")
    try:
        subprocess.check_call("makeblastdb -in %s -dbtype prot > /dev/null 2>&1" % f, shell=True)
    except:
        print("problem found in formatting annotation %s" % f)
    devnull = open('/dev/null', 'w')
    cmd = ["blastp",
           "-query", peptides,
           "-db", f,
           "-seg", my_seg,
           "-comp_based_stats", "F",
           "-num_threads", "1",
           "-evalue", "0.1",
           "-outfmt", "6",
           "-out", "%s_blast.out" % name]
    subprocess.call(cmd, stdout=devnull, stderr=devnull)
    devnull.close()

def _perform_workflow_blastn(data):
    tn = data[0]
    f = data[1]
    my_seg = data[2]
    peptides = data[3]
    algorithm = data[4]
    if ".fasta.new" in f:
        try:
            subprocess.check_call("makeblastdb -in %s -dbtype nucl > /dev/null 2>&1" % f, shell=True)
        except:
            print("problem found in formatting genome %s" % f)
    if ".fasta.new" in f:
        devnull = open('/dev/null', 'w')
        try:
            cmd = ["blastn",
                   "-task", algorithm,
                   "-query", peptides,
                   "-db", f,
                   "-dust", str(my_seg),
                   "-num_threads", "1",
                   "-evalue", "0.1",
                   "-outfmt", "6",
                   "-out", "%s_blast.out" % f]
            subprocess.call(cmd, stdout=devnull, stderr=devnull)
            devnull.close()
        except:
            print("The genome file %s was not processed" % f)

def blast_against_each_genome_blastn_dev(processors,algorithm,filter,peptides):
    """BLAST all peptides against each genome"""
    if "F" in filter:
        my_seg = "yes"
    else:
        my_seg = "no"
    curr_dir=os.getcwd()
    files = []
    for file in os.listdir(curr_dir):
        if file.endswith(".fasta.new"):
            files.append(file)
    files_and_temp_names = []
    for idx, f in enumerate(files):
        files_and_temp_names.append([str(idx), os.path.join(curr_dir,f), my_seg, peptides, algorithm])
    mp_shell(_perform_workflow_blastn, files_and_temp_names, processors)

def _perform_workflow_fdd(q, my_dict_o, data):
    tn = data[0]
    f = data[1]
    ref_scores = data[2]
    length = data[3]
    max_plog = data[4]
    min_hlog = data[5]
    clusters = data[6]
    processors = data[7]
    if "_blast.out" in f:
        genome_specific_dict = {}
        name = get_seq_name(f)
        reduced_name = name.replace(".fasta.new_blast.out","")
        genome_specific_dict["ID"] = reduced_name
        outfile = open("%s.counts.txt" % reduced_name, "w")
        try:
            with open(f) as infile:
                for line in infile:
                    fields = line.split()
                    if fields[0] not in ref_scores:
                        pass
                    elif float(fields[2]) >= int(min_hlog) and (float(fields[11]) / float(ref_scores.get(fields[0]))) >= float(length):
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
        for k,v in genome_specific_dict.items():
            if k == "ID":
                pass
            else:
                if k in clusters:
                    try:
                        new_dict.update({k:len(v)})
                    except:
                        new_dict.update({k:"0"})
        """This makes sure that every cluster is included"""
        od = collections.OrderedDict(sorted(new_dict.items()))
        ids = collections.OrderedDict({"ID":reduced_name})
        both = collections.OrderedDict(list(ids.items())+list(new_dict.items()))
        outfile.write(list(both.values())[0]+"\n")
        """This makes sure that the output is in order"""
        for cluster in clusters:
            if cluster in both:
                outfile.write(str(both.get(cluster))+"\n")
            else:
                outfile.write(str("0")+"\n")
        outfile.close()
        return genome_specific_dict

def find_dups_dev(ref_scores, length, max_plog, min_hlog, clusters, processors):
    from multiprocessing import Manager, Pool
    m = Manager()
    q = m.Queue()
    my_dict_o = m.dict()
    p = Pool(processors)
    curr_dir=os.getcwd()
    dup_dict = {}
    duplicate_file = open("duplicate_ids.txt", "w")
    genome_specific_list_of_lists = []
    files = os.listdir(curr_dir)
    files_and_temp_names = []
    for idx, f in enumerate(files):
        files_and_temp_names.append([str(idx), os.path.join(curr_dir, f), ref_scores, length, max_plog, min_hlog, clusters, processors])
    # Multiprocessing here (mp_shell for Ctrl+F)
    """How to test this function???"""
    for process in files_and_temp_names:
        p.apply(_perform_workflow_fdd, args=(q,my_dict_o,process))
    # Get rid of any duplicate values in queue
    unique = set()
    while q.empty() == False:
        unique.add(q.get())
    """This generates the list of all possible CDSs"""
    with open("dup_refs.txt", "a") as ref_file:
        ref_file.write("ID"+"\n")
        ref_file.write("\n".join(clusters)+"\n")
    ref_file.close()
    try:
        generate_dup_matrix()
        os.system("paste dup_refs.txt dup_values > dup_matrix.txt")
    except:
        print("problem generating duplicate matrix")
    """new way to report duplicates"""
    duplicate_IDs = []
    for line in open("dup_matrix.txt"):
        fields = line.split()
        if fields[0] == "ID":
            pass
        else:
            for field in fields[1:]:
                if float(field)>1:
                    if fields[0] in duplicate_IDs:
                        pass
                    else:
                        duplicate_IDs.append(fields[0])
    duplicate_file.write("\n".join(duplicate_IDs))
    duplicate_file.close()
    return duplicate_IDs

def _perform_workflow_nl(data):
     tn,f = data[0]
     clusters = data[1]
     names = data[2]
     table_list = data[3]
     name,values=make_table_test(f, "F", clusters)
     names.append(name)
     table_list.append(values)

def new_loop_dev(to_iterate, processors, clusters):
    from multiprocessing import Manager

    manager = Manager()
    names = manager.list()
    table_list = manager.list()

    files_and_temp_names_nl = []
    for file in to_iterate:
        files_and_temp_names_nl.append([file, clusters, names, table_list])

    mp_shell(_perform_workflow_nl, files_and_temp_names_nl, processors)

    names = list(names)
    table_list = list(table_list)

    return names,table_list

def make_table_test(infile, test, clusters):
    """make the BSR matrix table"""
    values = []
    names = []
    outdata = []
    name = get_seq_name(infile)
    reduced=[]
    """remove the junk at the end of the file"""
    reduced.append(name.replace('.fasta.new_blast.out.filtered.unique',''))
    names.append(reduced)
    my_dict={}
    with open(infile) as my_file:
        """make a dictionary of all clusters and values"""
        try:
            for line in my_file:
                fields=line.split()
                my_dict.update({fields[0]:fields[1]})
        except:
            raise TypeError("abnormal number of fields")
    """add in values, including any potentially missing ones"""
    for x in clusters:
        if x not in my_dict:
            my_dict.update({x:0})
    for x in reduced:
        values.append(x)
    """sort keys to get the same order between samples"""
    od = collections.OrderedDict(sorted(my_dict.items()))
    values += od.values()
    if "T" in test:
        return sorted(outdata)
    return names, values

def inverse_coding_regions(infile,ID):
    # Key = name of genome
    # Value = list of tuples, where each tuple is (start_range, stop_range)
    ranges = {}
    with open(infile) as f:
    # Get all ranges in input file
        for line in f:
            # Ignore metadata
            if line.startswith("#"):
                continue
            fields = line.split()
            try:
                name = fields[0]
                start = int(fields[3])
                stop = int(fields[4])
            # Not enough fields
            except IndexError as e:
                print("Error - unrecognized format in .gff file")
                print()
                print(e)
            # Trying to convert non-int to int
            except ValueError as e:
                print("Error - unrecognized format in .gff file")
                print()
                print(e)
            if name not in ranges:
                ranges[name] = [(start, stop)]
            else:
                ranges[name].append((start, stop))
    #ID is the name of the file, not the contig/chromosome
    outfile = open("%s.ranges" % ID, "w")
    # Sort the ranges based on the start range
    for name in ranges:
        ranges[name].sort(key=lambda tup: tup[0])
    for name in ranges:
        # Include range from 1 - (first start range) if first start range != 1
        if ranges[name][0][0] != 1:
            outfile.write("%s\t%s\t%s\n" % (name, str(1), str(ranges[name][0][0]-1)))
        for i in range(1, len(ranges[name])):
            previous_stop = ranges[name][i-1][1]
            current_start = ranges[name][i][0]
            # Ignore nested ranges
            if previous_stop+1 >= current_start:
                continue
            # Print (previous stop range + 1) - (current start range 1)
            outfile.write("%s\t%s\t%s\n" % (name, str(previous_stop+1), str(current_start-1)))
    outfile.close()

def parse_ranges_file(genome,ranges_file,name,test):
    """Make tuple of ranges file"""
    ranges_tuple = ()
    sequence = []
    if "/" in name:
        name_fields = name.split("/")
        reduced_name = name_fields[-1]
    else:
        reduced_name = name
    outfile = open("%s.intergenics.seqs" % reduced_name, "w")
    with open(ranges_file) as infile:
        for line in infile:
            newline = line.strip()
            fields = newline.split()
            ranges_tuple=((fields[0],fields[1],fields[2]),)+ranges_tuple
    with open(genome) as genome_file:
        for record in SeqIO.parse(genome_file,"fasta"):
            for range_tuple in ranges_tuple:
                if record.id == range_tuple[0]:
                    start = int(range_tuple[1])-1
                    end = int(range_tuple[2])
                    if len(str(record.seq[int(start):int(end)]))>50:
                        """This ignores regions shorter than 50 nucleotides, I think that these
                        should be renamed based on start and end of each range"""
                        outfile.write(">%s_%s_%s" % (range_tuple[0],start,end) +"\n")
                        outfile.write(str(record.seq[start:end])+"\n")
                        if "true" in test:
                            sequence.append(str(record.seq[start:end]))
    outfile.close()
    return sequence

def find_data_type(in_fasta):
    nt = []
    aa = []
    with open(in_fasta) as my_fasta:
        try:
            for record in SeqIO.parse(my_fasta, "fasta"):
                """None of these characters are IUPACs"""
                if "F" in record.seq or "L" in record.seq or "I" in record.seq or "P" in record.seq or "Q" in record.seq or "E" in record.seq:
                    aa.append("1")
                else:
                    pass
        except:
            print("file cannot be parsed, is it in FASTA format?")
            sys.exit()
    if len(aa)==0:
        data_type = "nt"
    elif len(aa)>0:
        data_type = "aa"
    else:
        print("problem determining data type..exiting")
        sys.exit()
    return data_type

#Putting the logging in here
import sys
import time

OUTSTREAM = sys.stdout
ERRSTREAM = sys.stderr

##
# Debug mode on?
DEBUG = False

def logPrint(msg, stream=None):
    if stream is None:
        stream = OUTSTREAM
    stream.write('LOG: %s - %s\n' % (timestamp(), removeRecursiveMsg(msg)))
    stream.flush()


def errorPrint(msg, stream=None):
    if stream is None:
        stream = ERRSTREAM

    stream.write('ERROR: %s - %s\n' % (timestamp(), removeRecursiveMsg(msg)))
    stream.flush()

def debugPrint(fmsg, stream=None):
    """In this case msg is a function, so the work is only done if debugging is one"""
    if DEBUG:
        if stream is None:
            stream = ERRSTREAM

        stream.write('DEBUG: %s - %s\n' % (timestamp(), removeRecursiveMsg(fmsg())))
        stream.flush()

def timestamp():
    return time.strftime('%Y/%m/%d %H:%M:%S')


def removeRecursiveMsg(msg):
    """
    This takes a message and if it starts with something that looks like
    a message generated with these tools it chops it off.  Useful if using
    one of these logging functions to print output from a program using
    the same logging functions
    """
    if msg.startswith('ERROR: ') or msg.startswith('DEBUG: ') or msg.startswith('LOG: '):
        return msg.split(' - ', 1)[1]
    else:
        return msg


##
# These version strip off the right white spaces so they can be used in printing output from a program
#errorPrintS = lambda l, *args, **kwargs : errorPrint(l.rstrip(), *args, **kwargs)

#logPrintS = lambda l, *args, **kwargs : logPrint(l.rstrip(), *args, **kwargs)

#debugPrintS = lambda f, *args, **kwargs : debugPrint(lambda : f().rstrip(), *args, **kwargs)
