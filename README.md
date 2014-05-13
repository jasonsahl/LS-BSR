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

###See manual for run directions###

  
