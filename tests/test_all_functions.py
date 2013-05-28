#!/usr/bin/env python

"""Calculate the BSR value for all predicted ORFs
in a set of genomes in fasta format.  V3 - replaced
transeq with BioPython.  V4 - changed to true BSR
test.  V5 - fixed bug in how BSR was calculated.
V6 - changed gene caller from Glimmer to Prodigal"""

import unittest
from ls_bsr.util import *
import os

class Test1(unittest.TestCase):
    def test(self):
        self.assertEqual(get_seq_name("/path/to/test.fasta"), "test.fasta")
class Test2(unittest.TestCase):
    def test(self):
        self.assertEqual(translate_consensus("test_consensus.fasta"), 'MTSFP')
class Test3(unittest.TestCase):
    def test(self):
        self.assertEqual(filter_seqs("test_length.pep"), [62])
class Test4(unittest.TestCase):
    def test(self):
        self.assertEqual(parse_self_blast("blast.test"), {'Cluster2': '60.6', 'Cluster0': '30.2', 'Cluster1': '40.5'})
class Test5(unittest.TestCase):
    def test(self):
        self.assertEqual(parse_blast_report("."), ['Cluster0', '30.2'])
class Test6(unittest.TestCase):
    def test(self):
        self.assertEqual(get_unique_lines("."), ['Cluster0\t30.2'])
class Test7(unittest.TestCase):
    def test(self):
        os.system("rm test_blast.out.filtered.filtered.unique")
        self.assertEqual(make_table(".", 2), ['30.2', '40.5', '60.6', '15.2', 0, '30.6'])
class Test8(unittest.TestCase):
    def test(self):
        self.assertEqual(divide_values( "test_bsr.matrix", {'Cluster2': '60.6', 'Cluster0': '30.2', 'Cluster1': '40.5'}),
                     [[1.0, 0.5033112582781457], [1.0, 0.0], [1.0, 0.504950495049505]])
class Test9(unittest.TestCase):
    def test(self):
        self.assertEqual(translate_genes("genes_test.fasta", "."), 'MNPHLTEHPPVGDIDALLQDTWLQVISLRQGVTCAEGEGQAFWQRCVADIERVHQALKDAGHSEQSCQHIRYAQCALLDETVKGRGVQDDAYFVWCHSPLQAHFFNTLDAGSQLYERMRAVLREPAPDRAVLTCFHRVLMLGFLGGYASPAASEREQLIDQLSVQVPAFSVAPSRGILASAASRNRLGIWLRYWPVRLGLAALMVALLWWGLDHWLSGLLATLLPEPV')
class Test10(unittest.TestCase):
    def test(self):
        self.assertEqual(rename_fasta_header("genes_test.fasta", "tmp.out"), ['>centroid_gi|22123922|ref|NC_004088.1|_3285', '>centroid_gi|22123922|ref|NC_004088.1|_1575'])
class Test11(unittest.TestCase):
    def test(self):
        self.assertEqual(autoIncrement(), 4)

def main():
    os.system("rm tmp.* test_blast.out.filtered sample*.filtered.unique.tmp.matrix ref.list names.txt genes.pep consensus.pep BSR_matrix_values.txt")

if __name__ == "__main__":
    unittest.main()
    main()

