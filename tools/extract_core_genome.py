#!/usr/bin/env python

"""If given core seqs and a
set of genomes, extract the core
and build an alignment ready
for phylogenetics"""

from __future__ import division
import optparse
import subprocess
import os
import sys
try:
    from Bio import SeqIO
    from Bio.Blast import NCBIXML
    from Bio.Seq import Seq
except:
    print "BioPython is not in your PATH, but needs to be"
    sys.exit()
import glob
import collections
try:
    from igs.utils import functional as func
    from igs.utils import logging
    from igs.threading import functional as p_func
except:
    "you need to add LS-BSR to your PYTHONPATH!"
    sys.exit()

"""errors for input start"""
def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print "directory of fastas cannot be found"
        sys.exit()

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print 'genes file cannot be opened'
        sys.exit()

def test_blast(option, opt_str, value, parser):
    if "tblastn" in value:
        setattr(parser.values, option.dest, value)
    elif "blastn" in value:
        setattr(parser.values, option.dest, value)
    else:
        print "Blast option not supported.  Only select from tblastn, blat, or blastn"
        sys.exit()
"""end, code starts now"""

def get_seq_name(fasta_in):
    name = os.path.basename(fasta_in)
    return name

def split_multifasta(infile):
    dir_out = os.getcwd()
    for record in SeqIO.parse(open(infile, "U"), "fasta"):
        f_out = os.path.join(dir_out,record.id+'.fasta')
        SeqIO.write([record],open(f_out,'w'),"fasta")

def combined_seqs(dir_path):
    num_genomes = []
    handle = open("combined.seqs", "w")
    for infile in glob.glob(os.path.join(dir_path, '*.fasta')):
        names = get_seq_name(infile)
        reduced = names.replace('.fasta','')
        print >> handle, ">"+str(reduced)
        for record in SeqIO.parse(open(infile), "fasta"):
            print >> handle, record.seq
    handle.close()
    for record in SeqIO.parse("combined.seqs", "fasta"):
        num_genomes.append(record.id)
    return len(num_genomes), num_genomes

def run_blast(infile, blast):
    names = get_seq_name(infile)
    reduced = names.replace('.fasta','')
    cmd = ["blastall",
            "-p", "%s" % blast,
            "-i", infile,
            "-d", "combined.seqs",
            "-o", "%s.blast.out" % reduced,
            "-m", "7",
            "-q", "-4",
            "-r", "5",
            "-a", "2",
            "-b", "2000",
            "-v", "2000",
            "-F", "F"]
    subprocess.check_call(cmd)
    return reduced

def parse_blast_xml_report(infile):
    """uses biopython to split the output
    from a blast file with xml output"""
    names = get_seq_name(infile)
    reduced = names.replace('.blast.out','')
    result_handle = open(infile, "U")
    blast_records = NCBIXML.parse(result_handle)
    blast_record = blast_records.next()    
    handle = open("%s.blast.parsed" % reduced, "w") 
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            test = Seq(hsp.sbjct)
            if int(hsp.query_start)<int(hsp.query_end):
                print >> handle, ">", alignment.title, test
            if int(hsp.query_start)>int(hsp.query_end):
                print >> handle, ">", alignment.title, test.reverse_complement()
    handle.close()
    result_handle.close()
    os.system("sort -u -k 3,3 %s.blast.parsed > %s.blast.unique" % (reduced, reduced))

def parsed_blast_to_seqs(infile):
    names = get_seq_name(infile)
    reduced = names.replace('.blast.unique','')
    outfile = open("%s.extracted.seqs" % reduced, "w")
    for line in open(infile, "U"):
        fields = line.split(" ")
        print >> outfile, fields[0] + fields[2], "\n", fields[3],
    outfile.close()

def check_and_align_seqs(infile, num_genomes):
    lengths = [ ]
    names = get_seq_name(infile)
    reduced = names.replace('.extracted.seqs','')
    for record in SeqIO.parse(infile, "fasta"):
        lengths.append(len(record.seq))
    lengths.sort(key=int)
    try:
        if (lengths[0]/lengths[-1]) > 0.75 and len(lengths) == num_genomes:
            os.system("muscle -in %s -out %s_aln.seqs > /dev/null 2>&1" % (infile,reduced))
        else:
            pass
    except:
        pass
    
def pull_seqs(names):
    curr_dir = os.getcwd()
    for infile in glob.glob(os.path.join(curr_dir, '*_aln.seqs')):
        for name in names:
            handle = open("%s_aln_final.seqs" % name, "a")
            for record in SeqIO.parse(open(infile), "fasta"):
                if name == record.id:
                    print >> handle, ">"+str(record.id)
                    print >> handle, record.seq
            handle.close()

def concatenate():
    curr_dir = os.getcwd()
    for infile in glob.glob(os.path.join(curr_dir, "*_aln_final.seqs")):
        names = get_seq_name(infile)
        reduced = names.replace("_aln_final.seqs", "")
        handle = open("%s.concat" % reduced, "w")
        print >> handle, "\n", ">"+str(reduced)
        for record in SeqIO.parse(open(infile), "fasta"):
            seqs = []
            seqs.append(record.seq)
            print >> handle, "".join([str(x) for x in seqs]),
        handle.close()

def fasta_to_tab(infile):
    my_file = open(infile, "rU")
    outfile = open("out.tab", "w")
    for record in SeqIO.parse(my_file, "fasta"):
        print >> outfile, record.id, record.seq
    my_file.close()
    outfile.close()

def tab_to_matrix(tab):
    reduced = [ ]
    out_matrix = open("tab_matrix", "w")
    for line in open(tab):
        tmp_list = []
        fields = line.split()
        tmp_list.append(fields[0])
        for nucs in fields[1]:
            tmp_list.append(nucs.upper())
        reduced.append(tmp_list)
    test=map(list, zip(*reduced))
    for x in test:
        print >> out_matrix, "\t".join(x)
    out_matrix.close()

def filter_alignment(tab):
    """currently untested, but needs to be"""
    outfile = open("tab.filtered", "w")
    infile = open(tab, "U")
    firstLine = infile.readline()
    print >> outfile, firstLine,
    for line in infile:
        valid_fields = []
        fields = line.split()
        if "-" in fields:
            pass
        else:
            print >> outfile, line,
    outfile.close()
    infile.close()

def file_to_fasta(matrix):
    """currently untested, but needs to be"""
    reduced = [ ]
    out_matrix = open("final_alignment.fasta", "w")
    for line in open(matrix, "U"):
        fields = line.strip().split()
        reduced.append(fields)
    test=map(list, zip(*reduced))
    for x in test:
        print >> out_matrix, ">"+str(x[0])
        print >> out_matrix, "".join(x[1:])
    out_matrix.close()
    
def remove_gaps(infile):
    fasta_to_tab(infile)
    tab_to_matrix("out.tab")
    filter_alignment("tab_matrix")
    file_to_fasta("tab.filtered")

def main(directory, genes, blast, processors, remove_gap, keep):
    dependencies = ['blastall','formatdb','muscle']
    for dependency in dependencies:
        ra = subprocess.call(['which', '%s' % dependency])
        if ra == 0:
            pass
        else:
            print "%s is not in your path, but needs to be!" % dependency
            sys.exit()
    start_dir = os.getcwd()
    ap=os.path.abspath("%s" % start_dir)
    dir_path=os.path.abspath("%s" % directory)
    try:
        os.makedirs('%s/to_extract_xxx' % ap)
        os.makedirs('%s/work_xxx' % ap)
    except:
        os.system("rm -rf %s/to_extract_xxx" % ap)
        os.system("rm -rf %s/work_xxx" % ap)
        os.makedirs('%s/to_extract_xxx' % ap)
        os.makedirs('%s/work_xxx' % ap)
    gene_path=os.path.abspath("%s" % genes)
    os.system("cp %s %s/to_extract_xxx/genes.fasta" % (gene_path,ap))
    os.chdir("%s/to_extract_xxx" % ap)
    split_multifasta("genes.fasta")
    os.system("rm genes.fasta")
    os.chdir("%s/work_xxx" % ap)
    """create combined file"""
    num_genomes, names = combined_seqs(dir_path)
    os.system("formatdb -i combined.seqs -p F")
    table_files = glob.glob(os.path.join("%s/to_extract_xxx" % ap, "*.fasta"))
    files_and_temp_names = [(str(idx), os.path.join("%s/to_extract_xxx" % ap, f))
                            for idx, f in enumerate(table_files)]
    def _perform_workflow(data):
        tn, f = data
        name = run_blast(f, blast)
        parse_blast_xml_report("%s.blast.out" % name)
        parsed_blast_to_seqs("%s.blast.unique" % name)
        check_and_align_seqs("%s.extracted.seqs" % name, num_genomes)
        os.system("rm %s.blast.out %s.blast.unique %s.extracted.seqs" % (name,name,name))
    set(p_func.pmap(_perform_workflow,
                    files_and_temp_names,
                    num_workers=processors))
    os.system("rm *.blast.out *.blast.unique *.extracted.seqs")
    pull_seqs(names)
    concatenate()
    os.system("cat *.concat > all.concat")
    os.system('sed "s/ //g" all.concat > tmp.concat')
    os.system("awk 'FNR>1' tmp.concat > all.concat")
    if remove_gap == "T":
        remove_gaps("all.concat")
        os.system("cp final_alignment.fasta %s" % ap)
    elif remove_gap == "F":
        os.system("cp all.concat %s/final_alignment.fasta" % ap)
    else:
        print "You have chosen an incorrect option for gap removal, choose from T or F"
        sys.exit()
    """finish up"""
    os.chdir("%s" % ap)
    if keep == "T":
        pass
    elif keep == "F":
        os.system("rm -rf %s/to_extract_xxx %s/work_xxx" % (ap,ap))
    else:
        print "Illegal keep value selected, not doing anything"
        pass
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-d", "--directory", dest="directory",
                      help="/path/to/fasta_directory [REQUIRED]",
                      type="string", action="callback", callback=test_dir)
    parser.add_option("-g", "--genes", dest="genes", action="callback", callback=test_file,
                      help="core genome genes to extract and align [REQUIRED]",
                      type="string")
    parser.add_option("-b", "--blast", dest="blast", action="callback", callback=test_blast,
                      help="blast method to use, blastn or tblastn, defaults to blastn",
                      type="string",default="blastn")
    parser.add_option("-p", "--parallel_workers", dest="processors",
                      help="How much work to do in parallel, defaults to 2",
                      default="2", type="int")
    parser.add_option("-r", "--remove_gaps", dest="remove_gap",
                      help="Should I remove gap columns? Defaults to T",
                      default="T", type="string", action="store")
    parser.add_option("-k", "--keep", dest="keep",
                      help="Keep temp files? Defaults to F",
                      default="F", type="string", action="store")
    options, args = parser.parse_args()
    mandatories = ["directory","genes"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.directory,options.genes,options.blast,options.processors,
         options.remove_gap,options.keep)
