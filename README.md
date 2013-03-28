LS-BSR is a method to compare all coding regions in a large set of genomes.
Each peptide is compared against it's nucleotide sequence in order to obtain
the maximum BLAST bit score.  Each peptide is then aligned against each genome
in order to find the query BLAST bit score.  The query dividied by the reference
provides one with the BSR, which can range from 0 to 1.

To run the program, the following dependencies are required:

1.  USearch (currently only works with v5)
2.  BioPython
3.  blastall (tested version is 2.2.25)
4.  glimmer3 (will be replaced with prodigal in future)

