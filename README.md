LS-BSR (Large Scale Blast Score Ratio) is released under the GPL version 3 license.  See "license.txt" for more information  

LS-BSR is a method to compare all coding regions in a large set of genomes.
Each peptide is compared against it's nucleotide sequence in order to obtain
the maximum BLAST bit score.  Each peptide is then aligned against each genome
in order to find the query BLAST bit score.  The query dividied by the reference
provides one with the BSR, which can range from 0 to 1; scores slightly higher
than 1.0 can be observed due to variable bit scores obtained by BLAST.  In my opinion,
they should be treated as 1.0.  Due to the "-C F" flag added recently, values > 1.00
have not been observed.  

contact: jasonsahl at gmail dot com  

Minimum requirements, see manual.md for version information  
1. Python >2.7 and <=3.5 (higher versions still work but tests may fail)
2. BioPython  
3. Prodigal - Required for de novo gene prediction only  
4. VSEARCH - Optional  
5. mmseqs2- Optional  
6. CD-HIT - Optional  
7. Blast+ - Optional  
8. Blat - Optional  
9. Diamond - Optional  

-To create an environment and run through conda:  
    `conda create -n ls_bsr python=3.9`  
    `conda activate ls_bsr`   
    `conda install -c bioconda blast vsearch cd-hit prodigal ucsc-blat diamond biopython mmseqs2`  
    #If you have problems with biopython, try: pip install Biopython  
    `git clone https://github.com/jasonsahl/LS-BSR.git`  
    `python setup.py install`  

-To test the install:  
    `python ls_bsr.py --version`  
    `python tests/test_all_functions.py`   

*See changelog.md for a list of changes  
*See manual.md for run directions  
