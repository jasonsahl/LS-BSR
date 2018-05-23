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
