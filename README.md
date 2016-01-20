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

##Minimum requirements, see manual for version information##
1. Python >2.7 and <3
2. BioPython
3. Prodigal - Required for de novo gene prediction only
4. VSEARCH - Optional
5. USEARCH - Optional
6. CD-HIT - Optional
7. Blast+ - Optional
8. Blat - Optional
9. NumPy - Optional

###Changes made on 1/5/2016###
1. A potential problem was observed when using USEARCH with the way the LS-BSR was clustering groups of clusters.
   This may have affected the overall composition of the pan-genome, but would not affect the associated BSR values
2. A bug was fixed in the way that VSEARCH and CD-HIT were receiving the input file for clustering. This should
   not have an effect on any results

###Significant changes made on 12/28/2015###
1. Genbank files now accepted as input, must end in ".gbk"
2. Your CDSs are only renamed if there is a conflicting sequence header. Otherwise,
   the original ID is used. If Prodigal prediction is selected, this will be the CDS ID that Prodigal assigns
3. CD-hit is now a supported clustering method
4. The "-u" and "-v" flags have been removed. The "-c" flag now is used to select the clustering method
   Choose from "usearch", "vsearch", or "cd-hit". These must now be in your $PATH variable as "usearch", "vsearch",
   or "cd-hit-est", respectively


###Significant changes made on 6/24/2015###

1. The use of blastall is deprecated. LS-BSR now uses blast+ instead.
2. VSEARCH can now be used instead of USEARCH. Notice that VSEARCH will likely be slower than USEARCH.

###See manual for run directions###
