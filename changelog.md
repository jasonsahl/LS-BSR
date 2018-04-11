
###Changes made on 4/11/2018###
1. LS-BSR is now Python3.5 compatible. Any Errors in testing shouldn't affect performance
2. LS-BSR now checks "-g" input to see if it is protein or nucleotide and catches errors accordingly
3. Duplicate searching can now be turned off with the "-z" flag
4. Changes to several of the external scripts for improved functionality

###Changes made on 10/23/2017###
1. Support for both blastp and diamond. These alignments are protein/protein and will be performed back against the gene predictions
2. I added the "-m" flag to Prodigal so that coding regions are not called through stretches of Ns. I also added the "-c" flag, which should prevent the creation of truncated coding regions that could inflate the pan-genome in some cases.

###Changes made on 8/30/2017###
1. intergenic sequences can now be compared using nucleotide aligners by using the "-y T" flag
2. The debug option ("-z") has been removed
3. A script has been included to convert LS-BSR output to Scoary (https://github.com/AdmiralenOla/Scoary) input

###Changes made on 2/9/2016###
1. BLASTN was changed to use BLASTN as the "task", instead of megablast. Using megablast can result in no alignments for more distant alignments, which would result in a BSR value of 0, even though significant alignments likely exist. This problem would only affect runs using BLASTN as the alignment method.
2. Because of the changes to BLASTN, the "penalty" and "reward" arguments have been removed and replaced with the defaults.


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

###Changes made on 12/13/2016###

1. Numpy is no longer required for any tool in LS-BSR
2. The unix paste function for the duplicate matrix creation has been replaced with a python function
3. All print statements have been replaced with print functions (still not entirely Python3 compatible, but on its way)
4. Added "-a", minimum peptide length. The previous default was set to "30", while new default is "33"
5. If nt sequences aren't in multiples of 3, they are trimmed back

###Changes made on 3/21/2017###

1. The "-y" flag has been removed. Now, the temporary folder is created in the current directory. If no prefix "-x" is provided,
   the temporary folder will be in $cwd/tmp. If you provide a prefix, the temporary directory will be "$CWD/prefix"
2. The multithreaded function has been changed for many of the functions. This shouldn't affect the function, but may result in less hanging.
3. The "paralogs.txt" file is no longer being produced. All duplicates are under the "duplicate_ids.txt"

