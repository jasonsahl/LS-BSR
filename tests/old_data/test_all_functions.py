#!/usr/bin/env python

"""Calculate the BSR value for all predicted ORFs
in a set of genomes in fasta format.  V3 - replaced
transeq with BioPython.  V4 - changed to true BSR
test.  V5 - fixed bug in how BSR was calculated.
V6 - changed gene caller from Glimmer to Prodigal"""

import unittest
from ls_bsr.util import *
import os
import tempfile
import shutil

curr_dir=os.getcwd()

class Test1(unittest.TestCase):
    def test(self):
        self.assertEqual(get_seq_name("/path/to/test.fasta"), "test.fasta")
    """tests the condition where you use a tilda instead of full path"""
    def test2(self):
        self.assertEqual(get_seq_name("~/test.fasta"), "test.fasta")
    """tests the case where no path is passed"""
    def test3(self):
        self.assertEqual(get_seq_name(""), "")
        
class Test2(unittest.TestCase):
    def test(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile")
        fp = open(fpath, "w")
        fp.write(">Cluster0\n")
        fp.write("ATGACGAGCTTTCCG")
        fp.close()
        self.assertEqual(translate_consensus(fpath), 'MTSFP')
        shutil.rmtree(tdir)
    def test2(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile")
        fp = open(fpath, "w")
        fp.write(">Cluster0\n")
        """tests the case of a difference seqeunce.  Also, tests
        wheter a stop codon will be recognized and excluded"""
        fp.write("ATGAATCACTACTAA")
        fp.close()
        self.assertEqual(translate_consensus(fpath), 'MNHY')
        shutil.rmtree(tdir)
    def test3(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile")
        fp = open(fpath, "w")
        fp.write(">Cluster1\n")
        """having an integer should make the script throw a typeerror"""
        fp.write("AT1CGAGCTTTCCG")
        fp.close()
        self.assertRaises(TypeError, translate_consensus, fpath)
        shutil.rmtree(tdir)
        os.system("rm tmp.pep")
        
class Test3(unittest.TestCase):
    def test(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile")
        fp = open(fpath, "w")
        fp.write(">Cluster0\n")
        """this peptide is 50 AA and should pass through filter"""
        fp.write("LHGRSCRAAFVTFGSTGYFGATAHEPARTTPTNARRRTTANRNACAAPDR\n")
        fp.write(">Cluster1\n")
        """this peptide is 49 AA and should get filtered out"""
        fp.write("LHGRSCRAAFVTFGSTGYFGATAHEPARTTPTNARRRTTANRNACAAPD\n")
        fp.write(">Cluster2\n")
        """empty lines won't throw an error, but will get filtered out"""
        fp.write(" ")
        fp.close()
        self.assertEqual(filter_seqs(fpath),[50])
        os.system("rm consensus.pep")
        shutil.rmtree(tdir)
    
class Test4(unittest.TestCase):
    def test(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile")
        fp = open(fpath, "w")
        fp.write("Cluster0	Cluster0	100.00	15	0	0	1	15	1	15	1e-07	30.2\n")
        fp.write("Cluster1	Cluster1	100.00	15	0	0	1	15	1	15	1e-07	40.5\n")
        fp.write("Cluster2	Cluster2	100.00	15	0	0	1	15	1	15	1e-07	60.6")
        fp.close()
        self.assertEqual(parse_self_blast(open(fpath,"U")), {'Cluster2': '60.6', 'Cluster0': '30.2', 'Cluster1': '40.5'})
        shutil.rmtree(tdir)
    def test2(self):
        """tests the condition where too few fields are observed in the blast report.
        should throw an error"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile")
        fp = open(fpath, "w")
        fp.write(">Cluster0\n")
        """this file has too few fields"""
        fp.write("Cluster0	Cluster0	100.00	15	0	0	1	15	1	15	1e-07")
        fp.close()
        self.assertRaises(TypeError, parse_self_blast, open(fpath, "U"))
        shutil.rmtree(tdir)
        
class Test5(unittest.TestCase):
    def test(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile_blast.out")
        os.chdir("%s" % tdir)
        fp = open(fpath, "w")
        fp.write("Cluster0	Cluster0	100.00	15	0	0	1	15	1	15	1e-07	30.2")
        fp.close()
        self.assertEqual(parse_blast_report(), ['Cluster0', '30.2'])
        os.chdir("%s" % curr_dir)
        shutil.rmtree(tdir)
    def test2(self):
        ndir = tempfile.mkdtemp(prefix="filetest_",)
        os.chdir("%s" % ndir)
        fpath = os.path.join(ndir,"output_blast.out")
        fp = open(fpath, "w")
        fp.write("Cluster0	Cluster0	100.00	15	0	0	1	15	1	15")
        fp.close()
        self.assertRaises(TypeError, parse_blast_report)
        os.chdir("%s" % curr_dir)
        shutil.rmtree(ndir)

class Test6(unittest.TestCase):
    def test(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered")
        os.chdir("%s" % tdir)
        fp = open(fpath, "w")
        fp.write("Cluster0	30.2\n")
        fp.write("Cluster0	15.3\n")
        fp.close()
        self.assertEqual(get_unique_lines(), ['Cluster0\t30.2\n'])
        os.chdir("%s" % curr_dir)
        shutil.rmtree(tdir)
    def test2(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered")
        os.chdir("%s" % tdir)
        fp = open(fpath, "w")
        """if file is empty, you can't get an error
        but you can get an empty set"""
        fp.write("")
        fp.close()
        self.assertEqual(get_unique_lines(), [])
        os.chdir("%s" % curr_dir)
        shutil.rmtree(tdir)
        
class Test7(unittest.TestCase):
    def test(self):
        """tests to make sure that list is being populated correctly.  The second file
        is missing a value and the list should be populated with a 0"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.fasta.new_blast.out.filtered.filtered.unique")
        fpath2 = os.path.join(tdir,"testfile2.fasta.new_blast.out.filtered.filtered.unique")
        os.chdir("%s" % tdir)
        fp = open(fpath, "w")
        fp.write("Cluster0	30.2\n")
        fp.write("Cluster1	40.5\n")
        fp.write("Cluster2	60.6")
        fp.close()
        fp2 = open(fpath2, "w")
        fp2.write("Cluster0	15.2\n")
        fp2.write("Cluster2	30.6")
        fp2.close()
        self.assertEqual(make_table(2), ['30.2', '40.5', '60.6', '15.2', 0, '30.6'])
        os.chdir("%s" % curr_dir)
        shutil.rmtree(tdir)
    def test2(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.fasta.new_blast.out.filtered.filtered.unique")
        fpath2 = os.path.join(tdir,"testfile2.fasta.new_blast.out.filtered.filtered.unique")
        fpath3 = os.path.join(tdir, "testfile3.fasta.new_blast.out.filtered.filtered.unique")
        os.chdir("%s" % tdir)
        fp = open(fpath, "w")
        fp.write("Cluster0	30.2\n")
        fp.write("Cluster1	40.5\n")
        fp.write("Cluster2	60.6")
        fp.close()
        fp2 = open(fpath2, "w")
        fp2.write("Cluster0	15.2\n")
        fp2.write("Cluster2	30.6")
        fp2.close()
        """tests the case where there is no value attached to an entry
        Test still works, even though an error is thrown"""
        #fp3 = open(fpath3, "w")
        #fp3.write("Cluster0")
        #fp3.close()
        #self.assertEqual(make_table(2), ['30.2', '40.5', '60.6', '15.2', 0, '30.6'])
        os.chdir("%s" % curr_dir)
        shutil.rmtree(tdir)
        """need to test if the file is completely empty"""
        
class Test8(unittest.TestCase):
    def test(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered")
        fp = open(fpath, "w")
        fp.write("	sample1	sample2\n")
        fp.write("Cluster0	30.2	15.2\n")
        fp.write("Cluster1	40.5	0\n")
        fp.write("Cluster2	60.6	30.6")
        fp.close()
        self.assertEqual(divide_values(fpath, {'Cluster2': '60.6', 'Cluster0': '30.2', 'Cluster1': '40.5'}),
                     [[1.0, 0.5033112582781457], [1.0, 0.0], [1.0, 0.504950495049505]])
        shutil.rmtree(tdir)
    def test2(self):
        """tests if a condition has a missing value"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered")
        fp = open(fpath, "w")
        fp.write("	sample1	sample2\n")
        fp.write("Cluster0	30.2	15.2\n")
        fp.write("Cluster1	40.5	0\n")
        fp.write("Cluster2	60.6")
        fp.close()
        self.assertRaises(TypeError, divide_values, {'Cluster2': '60.6', 'Cluster0': '30.2', 'Cluster1': '40.5'})
        os.system("rm BSR_matrix_values.txt")
        shutil.rmtree(tdir)
        
class Test9(unittest.TestCase):
    def test(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered")
        fp = open(fpath, "w")
        fp.write(">gi|22123922|ref|NC_004088.1|_3285\n")
        fp.write("ATGAATCCTCACCTAACCGAACACCCCCCAGTCGGGGATATTGACGCCCTGTTGCAGGACACCTGGCTACAGGTGATCAGCCTGCGTCAAGGGGTAACCTGTGCCGAGGGCGAAGGGCAGGCATTCTGGCAGCGCTGTGTGGCGGACATTGAACGTGTCCATCAGGCGCTGAAAGACGCCGGTCACAGCGAGCAGAGTTGCCAGCACATCCGATACGCCCAATGTGCACTGCTGGATGAG\n")
        fp.write(">gi|22123922|ref|NC_004088.1|_1575\n")
        fp.write("ATGAAGCTAAATATCAAAGTTAATTGTTCTTATATCTGTGAACCCATACGTAAGCAA")
        fp.close()
        """tests to see if the translation is correct, and if shorter sequences
        get filtered out"""
        self.assertEqual(translate_genes(fpath), 'MNPHLTEHPPVGDIDALLQDTWLQVISLRQGVTCAEGEGQAFWQRCVADIERVHQALKDAGHSEQSCQHIRYAQCALLDE')
        os.system("rm genes.pep")
        shutil.rmtree(tdir)
        
class Test10(unittest.TestCase):
    def test(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered")
        fp = open(fpath, "w")
        fp.write(">gi|22123922|ref|NC_004088.1|_3285\n")
        fp.write("ATGCGGGTTGGCCCGGGTTG\n")
        fp.write(">gi|22123922|ref|NC_004088.1|_1575\n")
        fp.write("MNPHLTEHPPVGDIDALLQDTWLQVISLRQGVT")
        fp.close()
        self.assertEqual(rename_fasta_header(fpath, "tmp.out"), ['>centroid_gi|22123922|ref|NC_004088.1|_3285', '>centroid_gi|22123922|ref|NC_004088.1|_1575'])
        os.system("rm tmp.out")
        shutil.rmtree(tdir)
         
class Test11(unittest.TestCase):
    def test(self):
        self.assertEqual(autoIncrement(), 4)

if __name__ == "__main__":
    unittest.main()
    main()

