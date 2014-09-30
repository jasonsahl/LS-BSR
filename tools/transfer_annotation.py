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
    
def main(peptides,consensus,processors):
    devnull = open("/dev/null", "w")
    pep_path=os.path.abspath("%s" % peptides)
    consensus_path=os.path.abspath("%s" % consensus)
    os.system("sed 's/ /_/g' %s > query.peptides.xyx" % pep_path)
    subprocess.check_call("formatdb -i query.peptides.xyx", shell=True)
    if consensus_path.endswith(".pep"):
        os.system("blastall -p blastp -i %s -d query.peptides.xyx -m 8 -o xyx.blast.out.xyx -a %s > /dev/null 2>&1" % (consensus_path,processors))
    elif consensus_path.endswith(".fasta"):
        os.system("blastall -p blastx -i %s -d query.peptides.xyx -m 8 -o xyx.blast.out.xyx -a %s > /dev/null 2>&1" % (consensus_path,processors))
    else:
        print "input genes file is of incorrect format, choose from fasta or pep"
        sys.exit()
    os.system("sort -u -k 1,1 xyx.blast.out.xyx > xyx.blast.unique.xyx")
    transfer_annotation(consensus, "xyx.blast.unique.xyx")
    os.system("rm query.peptides.xyx")

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
    options, args = parser.parse_args()
    
    mandatories = ["peptides", "consensus"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.peptides,options.consensus,options.processors)
    
