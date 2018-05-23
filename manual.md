## The Large Scale Blast Score Ratio (LS-BSR) pipeline

### Updated 4/10/2018

#### Citation:
Jason W. Sahl, J. Gregory Caporaso, David A. Rasko, Paul S. Keim (2014). The large-scale
blast score ratio (LS-BSR) pipeline: a method to rapidly compare genetic content between
bacterial genomes. PeerJ PrePrints 2:e220v1.

Jason W. Sahl, J. Gregory Caporaso, David A. Rasko, Paul S. Keim (2014). The large-scale
blast score ratio (LS-BSR) pipeline: a method to rapidly compare genetic content between
bacterial genomes. PeerJ 2 (e332).

#### Contact:
Please address queries, concerns, improvements to jasonsahl at gmail dot compare

#### What does it do?
The LS-BSR pipeline was designed as a way to quickly compare the genetic content between
a large number of bacterial genomes. LS-BSR can calculate several pan-genome statistics in a
population and the output can be easily visualized with a variety of third-party tools.
Additionally, LS-BSR can be used to query a set of genes and intergenic regions against a
large set of genomes to identify gene distribution and conservation. LS-BSR was developed to
be easy to run and interpret.

#### Recent changes:
1. I recently repeated some analyses to compare LS-BSR with Roary and the results are [here](https://github.com/jasonsahl/LS-BSR/wiki/Some-thoughts-on-comparing-Roary-to-LS-BSR)
2. LS-BSR now supports the use of Diamond and BLASTP for protein-protein comparisons.
One major difference by using this method is that the alignment is performed against the
gene predictions and not back against the genome.
3. Default parameters for Prodigal have been changed to prevent reading through large
stretches of Ns.
4. LS-BSR can now compare intergenic regions if the “-y” flag is set to “T” and if a nucleotide
aligner is chosen.
5. All functionality should now be compatible with Python 2.7 through 3.5. Seems to work with
3.6, although some additional testing is required.
6. Added “-z” toggle that can turn off duplicate searching, which can speed up very large
analyses.

#### Installation:

The code is kept [here](https://github.com/jasonsahl/LS-BSR.git)

-You can clone the repository to your own system with git:  
```git clone https://github.com/jasonsahl/LS-BSR.git```

-Enter the directory, then:  
```python setup.py install```

-if your install directory is /Users/jsahl/LS-BSR, run:  
```export PYTHONPATH=/Users/jsahl/LS-BSR:$PYTHONPATH```

-You can add this to your .bashrc or .profile  
-You can test your installation by running the tests:  
```python /Users/jsahl/LS-BSR/tests/test_all_functions.py```  
-If your installation is correct, all tests should pass  

#### Dependencies:  
1. Python >2.7 and <=3.5  
2. USEARCH (current tested version is 10.0.240): must be in your $PATH as “usearch” – At least
one clustering method must be chosen if a set of genes is not supplied. 32-bit version
should be sufficient for most applications, including the analysis of 1000 E. coli genomes.
Tested successfully with versions v7.0.959, 7.0.1090, 8.1.1861, 9.0.2132, and v10.0.240. Can
be freely obtained for academics/non-profits from (http://www.drive5.com/usearch/). Very
different results have been observed with V8 compared to earlier versions.  
3. VSEARCH (tested version is 1.1.3, but works with 2.0.4 and 2.5.0): must be in your $PATH as
“vsearch” – At least one clustering method must be chosen if a set of genes is not
supplied. Can be freely obtained at: [https://github.com/torognes/vsearch]. Currently, VSEARCH
does not work with protein sequences and you will see a warning if you try to combine
VSEARCH with BLASTP or DIAMOND.  
4. CD-HIT (tested version is 4.6): must be in your path as “cd-hit-est” for nucleotides and “cd-hit” for peptides - At least one clustering method must be chosen if a set of genes is not
supplied. Does not support clustering ID lower than 0.7 for nucleotides. Can be freely obtained
from: [https://github.com/weizhongli/cdhit]  
5. BioPython, must be in $PYTHONPATH environmental variable. Can be freely obtained from:
[http://biopython.org/wiki/Main_Page]  
6. Blast+ (tested version is 2.2.28 up to 2.2.50), must be in path as “blastn, tblastn, blastp,
makeblastdb’ – only required if you are using BLASTN, BLASTP, or TBLASTN, and not
BLAT. Blast+ can be obtained from: [ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.28/]. Weird and incorrect bit scores have been observed with 2.2.31. If BLASTP is used, the coding
regions are clustered, translated, and aligned against proteins predicted by Prodigal. This is a
fundamental difference compared to how BLASTN, BLAT, and TBLASTN work within LS-BSR.  
7. BLAT (tested version is v. 35x1), must be in path as ‘blat’ – only required if you use choose
blat for your alignment method. Can be obtained from:
[http://hgdownload.cse.ucsc.edu/admin/exe/]  
8. Prodigal (tested version is 2.60), must be in path as 'prodigal' - only required if a set of gene sequences is not supplied or if you want to use blastp or diamond. Can be obtained from:
[https://code.google.com/p/prodigal/]  
9. Diamond (tested version is v0.9.10.111), must be in path as ‘diamond’ – only required for
diamond protein comparisons. Can be obtained from: [https://github.com/bbuchfink/diamond]  

#### Command line options:
**-d DIRECTORY: --directory=DIRECTORY**: the directory to your fasta files, all must end in ".fasta". Can either be complete genomes or draft assemblies. Scaffolds are discouraged. Genbank files are supported and must end in “*.gbk” [REQUIRED]  
**-i ID**: de-replication clustering value, defaults to 0.9 (range from 0.0-1.0). Low values (<0.8) are not supported if CD-HIT is chosen as the clustering method  
**-f FILTER**: whether to use BLAST filtering, default is "F" or filter, turn off with "T". Turning this to “T” should speed up the analysis, but may throw out highly repetitive sequences.  
**-p PROCESSORS**: number of processors to use, defaults to 2.  
**-g GENES**: if you have a list of genes to screen, supply a nucleotide fasta file (.fasta) or a peptide file (.pep). Each gene sequence must be in frame, or questionable results will be obtained (only true for TBLASTN). If this flag is not invoked, then the de novo gene prediction method is invoked  
**-c CLUSTER_METHOD**: determines which clustering method to choose. You can choose from
“usearch”, “vsearch”, or “cd-hit”. These must be in your path as “usearch”, “vsearch”, “cd-hit-est”, or “cd-hit” to use.  
**-b BLAST**: which alignment method to use. Default is 'tblastn', can be changed to 'blastn', ‘blastp’, ‘diamond’, or ‘blat’. Can be used with either a list of supplied genes or with the de novo method. Tblastn, blastp, and diamond are not compatible with “-y T” flag set below.  
**-l LENGTH**: minimum BSR value to be called a duplicate, defaults to 0.7. The BSR of the "duplicate" divided by the reference bit score must be greater than this value to be called a
duplicate  
**-m MAX_PLOG**: maximum value to be called a remote paralog, defaults to 0.85. If the BSR value
is greater than this value, then it is considered to be a highly similar paralog  
**-n MIN_HLOG**: minimum BLAST ID to be called a highly similar paralog, defaults to 75. If the
BLAST ID is below this value, it is considered a remote paralog  
**-t F_PLOG**: filter ORFs with a remote paralog from BSR matrix? Default is F (do not filter), values, can be T (filter paralogs) or F  
**-k KEEP**: keep or remove temp files, choose from T or F, defaults to False (F), choose from T
or F  
**-s FILTER_PEPS**: filter out short peps < 50AA during TBLASTN? Defaults to True T), choose
from T or F  
**-e FILTER_SCAFFOLDS**: filter any contig that contains an N? Defaults to F, choose from T or
F  
**-x PREFIX**: prefix name for output files, defaults to time/date. If prefix is given, the temporary directory will be named after the prefix  
**-a min_pep_length**: after translating sequences, peptides of a length smaller than this value
will be discarded, defaults to 33 (integer)  
**-y intergenics**: Include intergenic regions greater than 50nts in the analysis? Regions at the
end of contigs will not be included. Choose from T or F, defaults to (F)    
**-z DUP_TOGGLE: Performs duplicate searching, which can take a while in large datasets.
Choose from T or F, defaults to “T”**  

#### Test data – give LS-BSR a whirl on small datasets  
Test data is present in the test_data directory. This data consists of:  
1. Genomes (4 E.coli genomes from 4 different pathogenic variants). Genomes are:  
* H10407 - enterotoxigenic E. coli (ETEC)  
* E2348/69 - enteropathogenic E. coli (EPEC)  
* O157:H7 sakai - shiga toxin E. coli (STEC)  
* SSON046 - Shigella sonnei  
2. Genes (5 different markers that delineate between the variants). These include:
* IpaH3 - Shigella invasion antigen. Mostly present in Shigella spp.  
* LT - heat-labile toxin. Only present in ETEC (not all)  
* ST2 - heat-stable toxin. Only present in ETEC (not all)  
* bfpB - bundle forming pilus. Only present on plasmid in EPEC  
* stx2a - shiga toxin. Present in STEC  

-You can test out the LS-BSR functionality in 4 different ways:    
1. Test the gene screen method with TBLASTN. Enter the test directory and run LS-BSR:  
```python /Users/jsahl/LS-BSR/ls_bsr.py -d genomes -g genes/ecoli_markers.fasta –x test```  
*the output (test_bsr_matrix.txt) should show how each gene is only present in the correct
pathovar  
2. Test the gene screen method with BLASTN. Enter the test directory and run LS-BSR:  
```python /Users/jsahl/LS-BSR/ls_bsr.py -d genomes -g genes/ecoli_markers.fasta -b blastn –x test```  
3. Test the de novo gene prediction method with USEARCH and TBLASTN for alignment:  
```python /Users/jsahl/LS-BSR/ls_bsr.py -d genomes –c usearch –x test```  
4. Test the diamond protein search method against annotations:  
```python ls_bsr.py -d genomes/ -b diamond -g genes.pep -x diamond```  

*Sample output for each method is shown in the test data directory:  
1. $prefix_bsr_matrix.txt: This is the 2x2 matrix of the BSR value for each CDS in each genome
queried  
2. $prefix_names.txt: The names of all of your genomes. This file can be helpful for running the
compare_BSR script described below  
3. $prefix_duplicate_ids.txt: A list of sequence IDs that are duplicated in at least one genome  
4. $prefix_consensus.fasta: A multi-fasta of all unique CDS sequences in the pan-genome  
5. $prefix_consensus.pep (optional): A multi-fasta of protein sequences if TBLASTN is selected  
6. $prefix_dup_matrix.txt: A 2x2 matrix showing how many copies of a CDS are present in each
genome, if conserved above a given threshold  
7. $prefix_run_parameters.txt: Documentation for the run that you performed  

#### Visualization of output:
1. The output of LS-BSR can be visualized in many different ways. One popular method to
visualize the output as a heatmap is through integration with R. Many beginners find R to be
intimidating, but this link by Kat Holt provides some excellent workflows on how to build
heatmaps, and correlate the output with phylogenies: [link](http://bacpathgenomics.wordpress.com/2012/05/25/displaying-data-associated-withphylogenetic-trees/)  
2. The Interactive Tree of Life (iTOL) project has an interface that is very user-friendly and can
directly take LS-BSR output and a phylogeny and create publication ready figures. iTOL can
be found [here](http://itol.embl.de/)  
3. MeV is designed as a way to visualize expression data, but can just as easily create
heatmaps from LS-BSR output. MeV can also be used to cluster LS-BSR data. MeV is
platform independent and can be found [here](http://www.tm4.org/)  
4. As part of the manuscript, we demonstrate how clustering the pan-genome can be used to
answer how genomes cluster by gene content. A new script shown below (BSR_to_cluster_dendrogram.py) can be used to quickly generate these cluster diagrams for further investigation.  

#### Post-matrix scripts:  
1. compare_BSR.py  
-what does it do? Looks for CDS differences between two user-defined populations.
Differences can be set by user-defined thresholds for presence and absence. The
“names.txt” file contains the names as they should be listed in your separate groups file
-what do you need for the script to run? Requirements include:  
• BSR matrix  
• Two new-line delimited group files, taken from “names.txt”  
• FASTA file of all CDS sequences  
-what does output look like? If there are unique sequences to either group, they will be
stored in the “groupX_unique_seqs.fasta” file. If there are no unique sequences, the
“groups_combined_header.txt” can be analyzed to look at the variable distribution of regions
between groups.  
```python compare_BSR.py -1 group1.txt -2 group2.txt –f $prefix_consensus.fasta –b $prefix_bsr_matrix.txt```  
2. filter_BSR_variome.py  
-what does it do? Filters out the conserved
