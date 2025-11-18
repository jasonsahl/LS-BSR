#!/usr/bin/env python

"""Calculate the BSR value for all predicted CDSs
in a set of genomes in fasta format.

written by Jason Sahl
contacted at jasonsahl@gmail.com
"""

from __future__ import print_function
import sys
import os
import optparse
from optparse import OptionParser
import subprocess
from subprocess import call
import errno
import types
from ls_bsr.util import *
import glob
import tempfile

def is_tool(name: str) -> bool:
    #from distutils.spawn import find_executable
    #return find_executable(name) is not None
    import shutil
    return shutil.which(name) is not None
    
def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s file cannot be opened' % option)
        sys.exit()

def test_filter(option, opt_str, value, parser):
    if "F" == value:
        setattr(parser.values, option.dest, value)
    elif "T" == value:
        setattr(parser.values, option.dest, value)
    else:
        print("option not supported.  Only select from T and F")
        sys.exit()

def test_cluster(option, opt_str, value, parser):
    if "mmseqs-lin" == value:
        setattr(parser.values, option.dest, value)
    elif "mmseqs" == value:
        setattr(parser.values, option.dest, value)
    elif "vsearch" == value:
        setattr(parser.values, option.dest, value)
    elif "cd-hit" == value:
        setattr(parser.values, option.dest, value)
    else:
        print("option not supported. Choose from vsearch, mmseqs, mmseqs-lin, cd-hit")
        sys.exit()

def test_blast(option, opt_str, value, parser):
    if "tblastn" == value:
        setattr(parser.values, option.dest, value)
    elif "blastn" == value:
        setattr(parser.values, option.dest, value)
    elif "blastn-short" == value:
        setattr(parser.values, option.dest, value)
    elif "blat" == value:
        setattr(parser.values, option.dest, value)
    elif "blastp" == value:
        setattr(parser.values, option.dest, value)
    elif "diamond" == value:
        setattr(parser.values, option.dest, value)
    else:
        print("Blast option not supported. Select from tblastn, blastn, blastn-short, blat, blastp, diamond")
        sys.exit()

def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print("directory of fastas cannot be found")
        sys.exit()

def test_id(option, opt_str, value, parser):
    if type(value) == int:
        sys.exit()
    elif type(value) == float:
        setattr(parser.values, option.dest, value)
    else:
        print("id value needs to be a float")
        sys.exit()

def test_fplog(option, opt_str, value, parser):
    if "F" in value:
        setattr(parser.values, option.dest, value)
    elif "T" in value:
        setattr(parser.values, option.dest, value)
    else:
        print("select from T or F for f_plog setting")
        sys.exit()

def main(directory,id,filter,processors,genes,cluster_method,blast,length,
         max_plog,min_hlog,f_plog,keep,filter_peps,filter_scaffolds,prefix,
         intergenics,min_len,dup_toggle):
    start_dir = os.getcwd()
    ap=os.path.abspath("%s" % start_dir)
    dir_path=os.path.abspath("%s" % directory)
    """Here's a check to make sure that there are no conflicting methods"""
    if "null" not in genes and "null" not in cluster_method:
        logPrint("Choose either genes or de novo clustering method, not both")
        sys.exit()
    """Test for use of intergenics with a protein alignment method"""
    if intergenics=="T" and blast=="tblastn":
        logPrint("Incompatible choices: if incorporating intergenics, choose a nucleotide alignment method")
        sys.exit()
    elif intergenics == "T" and blast=="blastp":
        logPrint("Incompatible choices: if incorporating intergenics, choose a nucleotide alignment method")
        sys.exit()
    elif intergenics == "T" and blast=="diamond":
        logPrint("Incompatible choices: if incorporating intergenics, choose a nucleotide alignment method")
        sys.exit()
    logPrint("Testing paths of dependencies")
    if blast=="blastn" or blast=="tblastn" or blast=="blastp":
        result = is_tool(blast)
        if result is False:
            print("blast isn't in your path, but needs to be!")
            sys.exit()
        else:
            print("citation: Altschul SF, Madden TL, Schaffer AA, Zhang J, Zhang Z, Miller W, and Lipman DJ. 1997. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Nucleic Acids Res 25:3389-3402")
    elif blast=="blat":
        result = is_tool("blat")
        if result is False:
            print("You have requested blat, but it is not in your PATH")
            sys.exit()
        else:
            print("citation: W.James Kent. 2002. BLAT - The BLAST-Like Alignment Tool.  Genome Research 12:656-664")
    elif blast=="diamond":
        result = is_tool("diamond")
        if result is False:
            print("You have requested DIAMOND, but it is not in your PATH (as diamond)")
            sys.exit()
        else:
            print("citation: Buchfink B, Xie C, Huson DH. 2015. Fast and sensitive protein alignment using DIAMOND. Nature methods, 12, 59-60.")
    elif blast=="blastn-short":
        print("citation: W.James Kent. 2002. BLAT - The BLAST-Like Alignment Tool.  Genome Research 12:656-664")
    if "NULL" in prefix:
        import datetime
        timestamp = datetime.datetime.now()
        tmp_rename = str(timestamp.year), str(timestamp.month), str(timestamp.day), str(timestamp.hour), str(timestamp.minute), str(timestamp.second)
        rename = "".join(tmp_rename)
        if os.path.exists("%s/%s" % (ap,rename)):
            print("old temp directory exists (%s/%s).  Delete and run again" % (ap,rename))
            sys.exit()
        else:
            os.makedirs("%s/%s" % (ap,rename))
            fastadir = ("%s/%s" % (ap,rename))
    else:
        if os.path.exists("%s/%s" % (ap,prefix)):
            print("old temp directory exists (%s/%s).  Delete and run again" % (ap,prefix))
            sys.exit()
        else:
            os.makedirs("%s/%s" % (ap,prefix))
            fastadir = "%s/%s" % (ap,prefix)
    samples = []
    for infile in glob.glob(os.path.join(dir_path, '*.fasta')):
        name=get_seq_name(infile)
        samples.append(name)
        os.link(infile,"%s/%s.new" % (fastadir,name))
    genbank_files = []
    for infile in glob.glob(os.path.join(dir_path, '*.gbk')):
        name=get_seq_name(infile)
        genbank_files.append("1")
    #New code to test if there are peptide files as the reference
    pep_refs = []
    for infile in glob.glob(os.path.join(dir_path, '*.pep')):
        name = get_seq_name(infile)
        pep_refs.append(name)
        pep_refs.append("1")
    if len(samples) == 0 and len(genbank_files) == 0 and len(pep_refs) == 0:
        print("no usable genome files found, exiting...")
        sys.exit()
    """This is the section on de novo clustering"""
    if "null" in genes:
        if "null" in cluster_method:
            print("Clustering method needed if genes aren't provided...exiting")
            sys.exit()
        else:
            pass
        result = is_tool("prodigal")
        if result is False:
            print("prodigal is not in your path, but needs to be!")
            sys.exit()
        else:
            print("citation: Hyatt D, Chen GL, Locascio PF, Land ML, Larimer FW, and Hauser LJ. 2010. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11:119")
        if "mmseqs" in cluster_method or "mmseqs-lin" in cluster_method:
            result = is_tool("mmseqs")
            if result is False:
                print("mmseqs is not in your path, but needs to be")
                sys.exit()
            else:
                print("citation: Steinegger and Soding. 2018. Clustering huge protein sequence sets in linear time. Nature Communications 9:2542")
        elif "cd-hit" in cluster_method:
            if blast == "blastp" or blast == "diamond" or blast == "tblastn":
                result = is_tool("cd-hit")
                if result is False:
                    print("cd-hit is not in your path, but needs to be!")
                    sys.exit()
                else:
                    print("citation: Li, W., Godzik, A. 2006. Cd-hit: a fast program for clustering and comparing large sets of protein or nuceltodie sequences. Bioinformatics 22(13):1658-1659")
        elif "vsearch" in cluster_method:
            if blast == "blastp" or blast == "diamond":
                print("vsearch not compatible with proteins, exiting...")
                sys.exit()
            else:
                result = is_tool("vsearch")
                if result is False:
                    print("vsearch is not in your path, but needs to be!")
                    sys.exit()
                else:
                    print("citation: Rognes, T., Flouri, T., Nichols, B., Qunice, C., Mahe, Frederic. 2016. VSEARCH: a versatile open source tool for metagenomics. PeerJ Preprints. DOI: https://doi.org/10.7287/peerj.preprints.2409v1")
        if len(samples) == 0:
            pass
        else:
            logPrint("predicting genes with Prodigal")
            """Only predict genes if there are FASTA files"""
            predict_genes(fastadir, processors, intergenics)
            logPrint("Prodigal done")
        """This function produces locus tags"""
        if len(genbank_files)>0:
            logPrint("Converting genbank files")
            os.chdir("%s" % fastadir)
            genbank_hits = process_genbank_files(dir_path)
        else:
            genbank_hits = []
        if genbank_files == None or len(genbank_files) == 0:
            if intergenics == "F":
                os.system("cat *genes.seqs > all_gene_seqs.out")
            elif intergenics == "T":
                os.system("cat *genes.seqs *intergenics.seqs > all_gene_seqs.out")
            else:
                pass
            if filter_scaffolds == "T":
                filter_scaffolds_fun("all_gene_seqs.out")
                os.system("mv tmp.out all_gene_seqs.out")
            else:
                pass
        else:
            """First combine all of the prodigal files into one file
            :First check that there will be prodigal annotations"""
            if len(samples)>0:
                if intergenics == "F":
                    os.system("cat *genes.seqs > all_gene_seqs.out.tmp")
                elif intergenics == "T":
                    os.system("cat *genes.seqs *intergenics.seqs > all_gene_seqs.out.tmp")
                if filter_scaffolds == "T":
                    filter_scaffolds_fun("all_gene_seqs.out.tmp")
                    os.system("mv tmp.out all_gene_seqs.out.tmp")
                else:
                    pass
            """This combines the locus tags with the Prodigal prediction.
            If there are no prodigal predictions, then no error is printed"""
            if len(samples)>0:
                os.system("cat *locus_tags.fasta all_gene_seqs.out.tmp > all_gene_seqs.out")
            else:
                os.system("cat *locus_tags.fasta > all_gene_seqs.out")
            if blast=="blastp" or blast=="diamond":
                """Need to convert the locus tags into peptides here"""
                translate_genes("all_gene_seqs.out","all_genes.pep",0)
                for infile in glob.glob(os.path.join(fastadir, "*locus_tags.fasta")):
                    base = os.path.basename(infile)
                    name = base.replace(".locus_tags.fasta","")
                    translate_genes(base,"%s.fasta.new_genes.pep" % name,0)
            else:
                for hit in genbank_hits:
                    reduced_hit = hit.replace(".gbk","")
                    """This is to ensure that genes are aligned back against the genome"""
                    SeqIO.convert("%s/%s" % (dir_path, hit), "genbank", "%s.fasta.new" % reduced_hit, "fasta")
        if "NULL" in cluster_method:
            print("Clustering chosen, but no method selected...exiting")
            sys.exit()
        elif "mmseqs" == cluster_method:
            logPrint("clustering with mmseqs at an ID of %s, using %s processors" % (id,processors))
            if blast == "blastp" or blast == "diamond":
                os.system("cat *new_genes.pep > all_gene_seqs.pep")
                run_mmseqs(id, processors, "all_gene_seqs.pep")
                os.system("mv mmseqs_rep_seq.fasta consensus.pep")
            else:
                run_mmseqs(id, processors, "all_gene_seqs.out")
                os.system("mv mmseqs_rep_seq.fasta consensus.fasta")
            logPrint("mmseqs clustering finished")
        elif "mmseqs-lin" == cluster_method:
            logPrint("clustering with mmseqs-linear at an ID of %s, using %s processors" % (id,processors))
            if blast == "blastp" or blast == "diamond":
                os.system("cat *new_genes.pep > all_gene_seqs.pep")
                run_mmseqs_lin(id, processors, "all_gene_seqs.pep")
                os.system("mv mmseqs_rep_seq.fasta consensus.pep")
            else:
                run_mmseqs_lin(id, processors, "all_gene_seqs.out")
                os.system("mv mmseqs_rep_seq.fasta consensus.fasta")
            logPrint("mmseqs-lin clustering finished")
        elif "vsearch" in cluster_method:
            logPrint("clustering with VSEARCH at an ID of %s, using %s processors" % (id,processors))
            run_vsearch(id, processors, "all_gene_seqs.out")
            os.system("mv vsearch.out consensus.fasta")
            logPrint("VSEARCH clustering finished")
        elif "cd-hit" in cluster_method:
            logPrint("clustering with cd-hit at an ID of %s, length percentage of %s, using %s processors" % (id,min_len,processors))
            if blast == "blastp" or blast == "diamond":
                os.system("cat *new_genes.pep > all_gene_seqs.pep")
                subprocess.check_call("cd-hit -i all_gene_seqs.pep -o consensus.pep -M 0 -T %s -c %s -s %s > cdhit.cluster 2>&1" % (processors,id,min_len), shell=True)
            else:
                subprocess.check_call("cd-hit-est -i all_gene_seqs.out -o consensus.fasta -M 0 -T %s -c %s -s %s > cdhit.cluster 2>&1" % (processors,id,min_len), shell=True)
        """need to check for dups here"""
        if os.path.exists("consensus.fasta"):
            dup_ids = test_duplicate_header_ids("consensus.fasta")
        elif os.path.exists("consensus.pep"):
            dup_ids = test_duplicate_header_ids("consensus.pep")
        else:
            print("clustering didn't work. Check input and try again")
            sys.exit()
        if dup_ids == "True":
            pass
        elif dup_ids == "False":
            print("duplicate headers identified, renaming..")
            try:
                rename_fasta_header("consensus.fasta", "tmp.txt")
                os.system("mv tmp.txt consensus.fasta")
            except:
                rename_fasta_header("consensus.pep", "tmp.txt")
                os.system("mv tmp.txt consensus.pep")
        else:
            pass
        if "tblastn" == blast:
            subprocess.check_call("makeblastdb -in consensus.fasta -dbtype nucl > /dev/null 2>&1", shell=True)
            translate_genes("consensus.fasta","consensus.pep",0)
            clusters = get_cluster_ids("consensus.pep")
            blast_against_self_tblastn("tblastn", "consensus.fasta", "consensus.pep", "tmp_blast.out", processors, filter)
        elif "blastn" == blast:
            subprocess.check_call("makeblastdb -in consensus.fasta -dbtype nucl > /dev/null 2>&1", shell=True)
            blast_against_self_blastn("blastn", "blastn", "consensus.fasta", "consensus.fasta", "tmp_blast.out", filter, processors)
            clusters = get_cluster_ids("consensus.fasta")
        elif "blat" == blast:
            blat_against_self("consensus.fasta", "consensus.fasta", "tmp_blast.out", processors)
            clusters = get_cluster_ids("consensus.fasta")
        elif "blastp" == blast:
            subprocess.check_call("makeblastdb -in consensus.pep -dbtype prot > /dev/null 2>&1", shell=True)
            blast_against_self_tblastn("blastp", "consensus.pep", "consensus.pep", "tmp_blast.out", processors, filter)
            clusters = get_cluster_ids("consensus.pep")
        elif "blastn-short" == blast:
            subprocess.check_call("makeblastdb -in consensus.fasta -dbtype nucl > /dev/null 2>&1", shell=True)
            blast_against_self_blastn("blastn", "blastn-short", "consensus.fasta", "consensus.fasta", "tmp_blast.out", filter, processors)
            clusters = get_cluster_ids("consensus.fasta")
        elif "diamond" == blast:
            subprocess.check_call("diamond makedb --in consensus.pep -d consensus > /dev/null 2>&1", shell=True)
            subprocess.check_call("diamond blastp -p 4 -d consensus -f 6 -q consensus.pep -o tmp_blast.out > /dev/null 2>&1", shell=True)
            clusters = get_cluster_ids("consensus.pep")
        subprocess.check_call("sort -u -k 1,1 tmp_blast.out > self_blast.out", shell=True)
        ref_scores=parse_self_blast("self_blast.out")
        os.system("cp tmp_blast.out ref.scores")
        subprocess.check_call("rm tmp_blast.out self_blast.out", shell=True)
        if "tblastn" == blast:
            logPrint("starting tblastn")
            blast_against_each_genome_tblastn_dev(processors, "consensus.pep", filter)
        elif "blastn" == blast:
            logPrint("starting blastn")
            blast_against_each_genome_blastn_dev(processors, "blastn", filter, "consensus.fasta")
        elif "blastn-short" == blast:
            logPrint("starting blastn-short")
            blast_against_each_genome_blastn_dev(processors, "blastn-short", filter, "consensus.fasta")
        elif "blat" == blast:
            logPrint("starting blat")
            blat_against_each_genome_dev("consensus.fasta",processors)
        elif "blastp" == blast:
            logPrint("starting blastp")
            blastp_against_each_annotation("consensus.pep",processors,filter)
        elif "diamond" == blast:
            logPrint("starting diamond")
            diamond_against_each_annotation("consensus.pep",processors)
        else:
            pass
    else:
        #########This section focuses on providing your own genes with -g############
        logPrint("Using pre-compiled set of predicted genes")
        files = glob.glob(os.path.join(dir_path,"*.fasta"))
        genbank_files = glob.glob(os.path.join(dir_path,"*.gbk"))
        pep_files = glob.glob(os.path.join(dir_path,"*.pep"))
        if len(genbank_files)>0:
            for hit in genbank_files:
                base = os.path.basename(hit)
                reduced_hit = base.replace(".gbk","")
                """This is to ensure that genes are aligned back against the genome"""
                os.chdir(fastadir)
                SeqIO.convert("%s/%s" % (dir_path, base), "genbank", "%s.fasta.new" % reduced_hit, "fasta")
                os.chdir(start_dir)
        if len(pep_refs)>0:
            for hit in pep_files:
                base = os.path.basename(hit)
                #This is creating problems. I'm commenting out to test the effect
                #reduced_hit = base.replace(".pep","")
                #os.link(hit,"%s/%s.new" % (fastadir,reduced_hit))
        if len(files)==0 and len(genbank_files)==0 and len(pep_refs) == 0:
            print("no usable reference genomes found!")
            sys.exit()
        gene_path = os.path.abspath("%s" % genes)
        data_type = find_data_type(gene_path)
        dup_ids = test_duplicate_header_ids(gene_path)
        #Test for duplicate headers and exit if found
        if dup_ids == "True":
            pass
        elif dup_ids == "False":
            print("duplicate headers identified, exiting..")
            sys.exit()
        clusters = get_cluster_ids(gene_path)
        os.chdir("%s" % fastadir)
        if gene_path.endswith(".pep"):
            if data_type == "aa":
                pass
            else:
                print("File is supposed to contain proteins, but needs to end with '.pep'..exiting")
                sys.exit()
            #Copies over the genes file into the working directory to prevent cross linking events
            os.system("cp %s %s/genes.pep" % (gene_path,fastadir))
            #Create the self alignment to get the reference scores
            if blast=="tblastn" or blast=="blastp":
                logPrint("using %s on peptides" % blast)
                try:
                    subprocess.check_call("makeblastdb -in genes.pep -dbtype prot > /dev/null 2>&1", shell=True)
                except:
                    logPrint("problem encountered formatting BLAST database")
                    sys.exit()
                blast_against_self_tblastn("blastp", "genes.pep", "genes.pep", "tmp_blast.out", processors, filter)
            elif blast=="diamond":
                logPrint("using %s on peptides" % blast)
                try:
                    subprocess.check_call("diamond makedb --in %s -d self > /dev/null 2>&1" % gene_path, shell=True)
                except:
                    logPrint("problem encountered formatting DIAMOND database")
                #Recently added masking flag
                subprocess.check_call("diamond blastp -p 4 -d self -f 6 -q %s -o tmp_blast.out --masking 0 > /dev/null 2>&1" % gene_path, shell=True)
            elif blast=="blat" or blast=="blastn" or blast=="blastn-short":
                print("Nucleotide aligner not compatible with protein sequences...exiting")
                sys.exit()
            #Aligning back against each genome
            if blast == "tblastn":
                logPrint("starting TBLASTN")
                blast_against_each_genome_tblastn_dev(processors,gene_path,filter)
            elif blast == "blastp":
                """I will need to first do gene prediction for each genome"""
                #First, check to see if the genomes are nt or pep
                if len(pep_refs)>0 and len(files)==0:
                    for infile in glob.glob(os.path.join(dir_path,'*.pep')):
                        name=get_seq_name(infile)
                        os.link(infile,"%s/%s.new" % (fastadir,name))
                else:
                    for infile in glob.glob(os.path.join(dir_path,'*.fasta')):
                        name=get_seq_name(infile)
                        os.link(infile,"%s/%s.new" % (fastadir,name))
                    #This indicates that there currently are no peptide files
                    if len(pep_files) == 0:
                        logPrint("Predicting genes with Prodigal")
                        predict_genes(fastadir, processors, intergenics)
                    else:
                        #TODO: add functionality here for proteomic screens
                        logPrint("performing peptide alignment against proteomes")
                logPrint("BlastP starting")
                #This script might need to be modified to fit with peptide "genomes"
                blastp_against_each_annotation(gene_path,processors,filter)
            elif blast == "diamond":
                for infile in glob.glob(os.path.join(dir_path, '*.pep')):
                #for infile in glob.glob(os.path.join(dir_path, '*.fasta')):
                    name=get_seq_name(infile)
                    os.link(infile,"%s/%s.new" % (fastadir,name))
                if len(pep_files)==0:
                    logPrint("Predicting genes with Prodigal")
                    predict_genes(fastadir, processors, intergenics)
                    logPrint("Diamond starting")
                    diamond_against_each_annotation(gene_path,processors)
                else:
                    logPrint("performing peptide alignment against proteomes")
                    #If the file ends in ".pep", then the name should be "genome.pep.new. currently , thare are two files with a "new" suffix
                    diamond_against_each_proteome(gene_path,processors,"%s/%s.new" % (fastadir,name))
        elif gene_path.endswith(".fasta"):
            if data_type == "nt":
                pass
            else:
                print("File is supposed to contain nucleotides, but doesn't look correct..exiting ")
                sys.exit()
            os.system("cp %s %s" % (gene_path,fastadir))
            if blast == "diamond" or blast == "blastp":
                print("protein alignment not compatible with nucleotide input..exiting")
                sys.exit()
            if "tblastn" == blast:
                logPrint("using tblastn")
                translate_genes(gene_path,"genes.pep",0)
                try:
                    subprocess.check_call("makeblastdb -in %s -dbtype nucl > /dev/null 2>&1" % gene_path, shell=True)
                except:
                    logPrint("problem encountered with BLAST database")
                    sys.exit()
                blast_against_self_tblastn("tblastn", gene_path, "genes.pep", "tmp_blast.out", processors, filter)
                logPrint("starting BLAST")
                blast_against_each_genome_tblastn_dev(processors, "genes.pep", filter)
                os.system("cp genes.pep %s" % start_dir)
            elif "blastn" == blast:
                logPrint("using blastn")
                try:
                    subprocess.check_call("makeblastdb -in %s -dbtype nucl > /dev/null 2>&1" % gene_path, shell=True)
                except:
                    logPrint("Database not formatted correctly...exiting")
                    sys.exit()
                try:
                    blast_against_self_blastn("blastn", "blastn", gene_path, gene_path, "tmp_blast.out", filter, processors)
                except:
                    print("problem with blastn, exiting")
                    sys.exit()
                logPrint("starting BLAST")
                try:
                    blast_against_each_genome_blastn_dev(processors,"blastn",filter,gene_path)
                except:
                    print("problem with blastn, exiting")
                    sys.exit()
            elif "blastn-short" == blast:
                logPrint("using blastn-short")
                try:
                    subprocess.check_call("makeblastdb -in %s -dbtype nucl > /dev/null 2>&1" % gene_path, shell=True)
                except:
                    logPrint("Database not formatted correctly...exiting")
                    sys.exit()
                try:
                    blast_against_self_blastn("blastn", "blastn-short", gene_path, gene_path, "tmp_blast.out", filter, processors)
                except:
                    print("problem with blastn, exiting")
                    sys.exit()
                logPrint("starting BLAST")
                try:
                    blast_against_each_genome_blastn_dev(processors,"blastn-short",filter,gene_path)
                except:
                    print("problem with blastn-short, exiting")
                    sys.exit()
            elif "blat" == blast:
                logPrint("using blat")
                blat_against_self(gene_path, gene_path, "tmp_blast.out", processors)
                logPrint("starting BLAT")
                blat_against_each_genome_dev(gene_path,processors)
            else:
                pass
        else:
            print("input file format not supported")
            sys.exit()
        subprocess.check_call("sort -u -k 1,1 tmp_blast.out > self_blast.out", shell=True)
        os.system("cp self_blast.out ref.scores")
        ref_scores=parse_self_blast("self_blast.out")
        subprocess.check_call("rm tmp_blast.out self_blast.out", shell=True)
        """testing block complete"""
    if blast=="blat":
        logPrint("BLAT complete")
    elif blast=="diamond":
        logPrint("Diamond complete")
    else:
        logPrint("BLAST done")
    if dup_toggle == "T" or "T" in f_plog:
        logPrint("Finding duplicates")
        find_dups_dev(ref_scores, length, max_plog, min_hlog, clusters, processors)
        logPrint("Finding duplicates complete")
    else:
        logPrint("Duplicate searching turned off")
    parse_blast_report_dev("false",4)
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
    logPrint("starting matrix building")
    new_names,new_table = new_loop_dev(files_and_temp_names, processors, clusters)
    new_table_list = table_list+new_table
    open("ref.list", "a").write("\n")
    for x in nr_sorted:
        open("ref.list", "a").write("%s\n" % x)
    names_out = open("names.txt", "w")
    names_redux = [val for subl in new_names for val in subl]
    #Testing this method
    for x in names_redux:
        names_out.write("".join(x)+"\n")
    names_out.close()
    create_bsr_matrix_dev(new_table_list)
    divide_values("bsr_matrix",ref_scores)
    subprocess.check_call("paste ref.list BSR_matrix_values.txt > %s/bsr_matrix_values.txt" % start_dir, shell=True)
    try:
        if dup_toggle == "T":
            subprocess.check_call("cp dup_matrix.txt names.txt consensus.pep duplicate_ids.txt consensus.fasta %s" % ap, shell=True, stderr=open(os.devnull, 'w'))
        else:
            subprocess.check_call("cp names.txt consensus.pep consensus.fasta %s" % ap, shell=True, stderr=open(os.devnull, 'w'))
    except:
        pass
    if "T" in f_plog:
        logPrint("filtering duplicates")
        num_filtered = filter_paralogs("%s/bsr_matrix_values.txt" % start_dir, "duplicate_ids.txt")
        logPrint("%s duplicates filtered" % str(num_filtered))
        if "NULL" in prefix:
            os.system("cp bsr_matrix_values_filtered.txt %s/%s_paralogs_filtered_bsr_matrix_values.txt" % (start_dir,"".join(rename)))
        else:
            os.system("cp bsr_matrix_values_filtered.txt %s/%s_paralogs_filtered_bsr_matrix_values.txt" % (start_dir,prefix))
    os.chdir("%s" % ap)
    logPrint("matrix built")
    if "NULL" in prefix:
        if dup_toggle == "T":
            os.system("mv dup_matrix.txt %s_dup_matrix.txt" % "".join(rename))
            os.system("mv duplicate_ids.txt %s_duplicate_ids.txt" % "".join(rename))
        else:
            pass
        os.system("mv names.txt %s_names.txt" % "".join(rename))
        os.system("mv bsr_matrix_values.txt %s_bsr_matrix.txt" % "".join(rename))
        if os.path.isfile("consensus.fasta"):
            os.system("mv consensus.fasta %s_consensus.fasta" % "".join(rename))
        if os.path.isfile("consensus.pep"):
            os.system("mv consensus.pep %s_consensus.pep" % "".join(rename))
    else:
        if dup_toggle == "T":
            os.system("mv dup_matrix.txt %s_dup_matrix.txt" % prefix)
            os.system("mv duplicate_ids.txt %s_duplicate_ids.txt" % prefix)
        else:
            pass
        os.system("mv names.txt %s_names.txt" % prefix)
        os.system("mv bsr_matrix_values.txt %s_bsr_matrix.txt" % prefix)
        if os.path.isfile("consensus.fasta"):
            os.system("mv consensus.fasta %s_consensus.fasta" % prefix)
        if os.path.isfile("consensus.pep"):
            os.system("mv consensus.pep %s_consensus.pep" % prefix)
    if "NULL" in prefix:
        outfile = open("%s_run_parameters.txt" % "".join(rename), "w")
    else:
        outfile = open("%s_run_parameters.txt" % prefix, "w")
    outfile.write("-d %s \\\n" % directory)
    outfile.write("-i %s \\\n" % id)
    outfile.write("-f %s \\\n" % filter)
    outfile.write("-p %s \\\n" % processors)
    outfile.write("-g %s \\\n" % genes)
    outfile.write("-c %s \\\n" % cluster_method)
    outfile.write("-b %s \\\n" % blast)
    outfile.write("-l %s \\\n" % length)
    outfile.write("-m %s \\\n" % max_plog)
    outfile.write("-n %s \\\n" % min_hlog)
    outfile.write("-t %s \\\n" % f_plog)
    outfile.write("-k %s \\\n" % keep)
    outfile.write("-s %s \\\n" % filter_peps)
    outfile.write("-e %s \\\n" % filter_scaffolds)
    outfile.write("-x %s \\\n" % prefix)
    outfile.write("-y %s \\\n" % intergenics)
    outfile.write("-ml %s \\\n" % min_len)
    outfile.write("-z %s\n" % dup_toggle)
    outfile.write("temp data stored here if kept: %s" % fastadir)
    outfile.close()
    logPrint("all Done")
    if "T" == keep:
        pass
    else:
        os.system("rm -rf %s" % fastadir)
    os.chdir("%s" % ap)

if __name__ == "__main__":
    parser = OptionParser(usage="usage: %prog [options]",version="%prog 1.8")
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
                      help="Clustering method to use: choose from mmseqs, mmseqs-lin, vsearch, cd-hit",
                      type="string", default="null")
    parser.add_option("-b", "--blast", dest="blast", action="callback", callback=test_blast,
                      help="use tblastn, blastn, blastp, blastn-short, diamond, or blat (nucleotide search only), default is tblastn",
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
                      help="remove short peptides, smaller than 50AA?  Defaults to F",
                      default="T", callback=test_filter, type="string")
    parser.add_option("-e", "--filter_scaffolds", dest="filter_scaffolds", action="callback",
                      help="Filter any contig that contains an N? Defaults to F",
                      default="F", callback=test_filter, type="string")
    parser.add_option("-x", "--prefix", dest="prefix", action="store",
                      help="prefix for naming output files, defaults to time/date",
                      default="NULL", type="string")
    parser.add_option("-y", "--intergenics", dest="intergenics", action="callback",
                      help="Incoporate intergenic regions? T or F; Defaults to F",
                      default="F", type="string", callback=test_filter)
    parser.add_option("--ml", "--min_len", dest="min_len", action="store",
                      help="value for 's' flag in cd-hit, defaults to '-i' flag in LS-BSR (float)",
                      default="0.9", type="float")
    parser.add_option("-z", "--dup_toggle", dest="dup_toggle", action="callback",
                      help="Perform duplicate searching? T or F; Defaults to F",
                      default="F", type="string", callback=test_filter)
    options, args = parser.parse_args()

    mandatories = ["directory"]
    for m in mandatories:
        if not getattr(options, m, None):
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)

    main(options.directory,options.id,options.filter,options.processors,options.genes,options.cluster_method,options.blast,
         options.length,options.max_plog,options.min_hlog,options.f_plog,options.keep,options.filter_peps,
         options.filter_scaffolds,options.prefix,options.intergenics,options.min_len,options.dup_toggle)
