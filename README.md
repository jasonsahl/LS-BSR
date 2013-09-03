LS-BSR (Large Scale Blast Score Ratio) is released under the MIT license.  See "license.txt" for more information

LS-BSR is a method to compare all coding regions in a large set of genomes.
Each peptide is compared against it's nucleotide sequence in order to obtain
the maximum BLAST bit score.  Each peptide is then aligned against each genome
in order to find the query BLAST bit score.  The query dividied by the reference
provides one with the BSR, which can range from 0 to 1; scores slightly higher
than 1.0 can be observed due to variable bit scores obtained by BLAST.  In my opinion,
they should be treated as 1.0.

contact: jasonsahl at gmail dot com

To run the program, the following dependencies are required:

1.  USearch (tested version is 6.0.307), path is passed as command-line option - only required
    if a set of gene sequences is not supplied.  32-bit version should be sufficient.  Tested successfully
    with version v7.0.959.
2.  BioPython, must be in PythonPath
3.  blastall (tested version is 2.2.25), must be in path as 'blastall'.  Known issues
    seen with v2.2.26
4.  Prodigal (tested version is 2.60), must be in path as 'prodigal' - only required
    if a set of gene sequences is not supplied
5.  Numpy, must be in PythonPath.  Numpy is only required for the compare matrices tool.
	If you don't want to install numpy, comment out this line in ls_bsr/util.py
	
	import numpy as np


##To install:

python setup.py install

-if your install directory is /Users/jsahl/LS-BSR, then run:

export PYTHONPATH=/Users/jsahl/LS-BSR:$PYTHONPATH

You can test your installation by running the tests:

python /Users/jsahl/LS-BSR/tests/test_all_functions.py

-If your installation is correct, all 54 tests should pass

Command line options include:

-d DIRECTORY, --directory=DIRECTORY : the directory to your fasta files, all must end in 
".fasta"

-i ID, de-replication clustering value for USEARCH, defaults to 0.9 (range from 0.0-1.0)

-f FILTER, whether to use BLAST filtering, default is "F" or filter, turn off with "T"

-p PROCESSORS, number of processors to use, defaults to 2

-g GENES, if you have a list of genes to screen, supply a nucleotide fasta file. Each gene
sequence must be in frame, or questionable results will be obtained.  If this flag is not envoked,
 then the de novo gene prediction method is applied
 
-b BLAST, which blast method to use for a supplied set of genes.  Default is 'tblastn',
  can be changed to 'blastn'.  Only is functional when -g flag is invoked
  
-q PENALTY, blast mismatch penalty, default is -4, only works with blastn and -g flag.
   Optimized to return longer matches.  Only certain q/r ratios are allowed.  See BLAST
   documentation for more details.

-r REWARD, blast reward value, default is 5, only works with blastn and -g flag.
   Optimized to return longer matches.  Only certain q/r ratios are allowed.  See BLAST
   documentation for more details.
   
-l LENGTH, minimum BSR value to be called a duplicate, defaults to 0.7.  The BSR of the
   "duplicate" divided by the reference bit score must be greater than this value to be
   called a duplicate
   
-m MAX_PLOG, maximum value to be called a paralog, defaults to 0.85.  If the BSR value
   is greater than this value, then it is considered to be an ortholog
  
-n MIN_HLOG, minimum BLAST ID to be called a homolog, defaults to 75.  If the BLAST
   ID is below this value, it is considered a remote homolog
   
-t F_PLOG, filter ORFs with a paralog from BSR matrix? Default is F (do not filter), 
   values can be T (filter paralogs) or F

   
Test data is present in the test_data directory.  This data consists of:

1.  Genomes (4 E.coli genomes from 4 different pathogenic variants).  Genomes are:

-H10407 - enterotoxigenic E. coli (ETEC) 
-E2348/69 - enteropathogenic E. coli (EPEC) 
-O157:H7 sakai - shiga toxin E. coli (STEC) 
-SSON046 - Shigella sonnei 

2.  Genes (5 different markers that deliniate between the variants).  These include:

IpaH3 - Shigella invasion antigen.  Mostly present in Shigella spp.
LT - heat-labile toxin.  Only present in ETEC
ST2 - heat-stable toxin.  Only present in ETEC
bfpB - bundle forming pilus.  Only present on plasmid in EPEC
stx2a - shiga toxin.  Present in STEC

You can test out the LS-BSR functionality in 3 different ways:

1.  Test the gene screen method with tblastn:

-enter test_data directory

-run LS-BSR

python /Users/jsahl/LS-BSR/ls_bsr.py -d genomes -g genes/ecoli_markers.fasta

-the output should show how each gene is only present in the correct pathovar

2. Test the gene screen method with blastn:

-enter test_data directory

-run LS-BSR

python /Users/jsahl/LS-BSR/ls_bsr.py -d genomes -g genes/ecoli_markers.fasta -b blastn

3.  Test the de novo gene prediction method:

-enter test_data directory

-run LS-BSR

python /Users/jsahl/LS-BSR/ls_bsr.py -d genomes -u /usr/local/bin/usearch6

-To inspect the output, you can look up the following entries in the BSR matrix.  They
should correspond with the results obtained with the gene screen methods:

IpaH3 -> centroid_1724,
LT -> centroid_11953,
ST2 -> centroid_19265,
bfpB -> centroid_1922,
stx2a -> centroid_7471

Sample output for each method is shown in the test_data directory

##Tools
1.  compare_BSR_matrices.py - Find CDSs conserved in one group and missing from another.
    Group names must be in new-line delimited text files.  The default is to identify CDSs
    with a BSR >=0.8 in one group and <0.4 in the second group.  These values can be
    changed with command line parameters
2.  filter_column_BSR.py - Filter out a sample from your BSR matrix.  This is helpful
    if you have a problematic genome
3.  pan_genome_stats.py - Calculates core and accessory genome statistics from user-defined
    thresholds
4.  filter_BSR_variome.py - Filters out the "core" genome based on a user-defined threshold.
    Sometimes we are only interested in the variable part of the pan-genome     
                   

  
