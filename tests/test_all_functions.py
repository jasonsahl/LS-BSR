"""test functions for LS-BSR tools"""

import unittest
from ls_bsr.util import *
import os
import tempfile
import shutil

curr_dir=os.getcwd()

class Test1(unittest.TestCase):
    def test_get_seq_name_basic_function(self):
        self.assertEqual(get_seq_name("/path/to/test.fasta"), "test.fasta")
    """tests the condition where you use a tilda instead of full path"""
    def test_get_seq_name_tilda(self):
        self.assertEqual(get_seq_name("~/test.fasta"), "test.fasta")
    """tests the case where no path is passed"""
    def test_get_seq_name_empty(self):
        self.assertEqual(get_seq_name(""), "")
    """tests the case where something weird is passed"""
    def test_get_seq_name_wrong_slash(self):
        self.assertEqual(get_seq_name("\wrong\way"), "\\wrong\\way")

class Test3(unittest.TestCase):
    def test_filter_seqs_length_filter(self):
        """Tests several conditions in one shot"""
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
        self.assertEqual(filter_seqs(fpath,"tmp.pep"),[50])
        os.system("rm tmp.pep")
        shutil.rmtree(tdir)
    def test_filter_seqs_non_nucleotide(self):
        """tests condition where you have non nucleotide
        characters.  Won't throw an error, but will give
        you an empty set"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile")
        fp = open(fpath, "w")
        fp.write(">Cluster0\n")
        fp.write("12344423432343")
        fp.close()
        self.assertEqual(filter_seqs(fpath,"tmp.pep"),[])
        os.system("rm tmp.pep")
        shutil.rmtree(tdir)

class Test4(unittest.TestCase):
    def test_parse_self_blast_basic_function(self):
        """Tests the basic functionality of the parse_self_blast
        function"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile")
        fp = open(fpath, "w")
        fp.write("Cluster0\tCluster0\t100.00\t15\t0\t0\t1\t15\t1\t15\t1e-07	\t30.2\n")
        fp.write("Cluster1	Cluster1	100.00	15	0	0	1	15	1	15	1e-07	40.5\n")
        fp.write("Cluster2	Cluster2	100.00	15	0	0	1	15	1	15	1e-07	60.6")
        fp.close()
        self.assertEqual(parse_self_blast(fpath), {'Cluster2': '60.6', 'Cluster0': '30.2', 'Cluster1': '40.5'})
        shutil.rmtree(tdir)
    def test_parse_self_blast_missing_fields(self):
        """tests the condition where too few fields are observed in the blast report.
        should throw an error"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile")
        fp = open(fpath, "w")
        """this file has too few fields"""
        fp.write("Cluster0	Cluster0	100.00	15	0	0	1	15	1	15	1e-07")
        fp.close()
        self.assertRaises(TypeError, parse_self_blast, fpath)
        shutil.rmtree(tdir)
    def test_parse_self_blast_empty_file(self):
        """tests the condition where the file is empty
        should create an empty dictionary"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"empty_file")
        fp = open(fpath, "w")
        fp.write("")
        fp.close()
        self.assertEqual(parse_self_blast(fpath),{})
        shutil.rmtree(tdir)

class Test8(unittest.TestCase):
    def test_divide_values_basic_function(self):
        """tests basic functionality of the divide_values function"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered")
        fp = open(fpath, "w")
        fp.write("\tsample1\tsample2\n")
        fp.write("Cluster0\t30.2\t15.2\n")
        fp.write("Cluster1\t40.5\t0\n")
        fp.write("Cluster2\t60.6\t30.6")
        fp.close()
        self.assertEqual(divide_values(fpath, {'Cluster2': '60.6', 'Cluster0': '30.2', 'Cluster1': '40.5'}),
                     [[1.0, 0.5033112582781457], [1.0, 0.0], [1.0, 0.504950495049505]])
        shutil.rmtree(tdir)
    def test_divide_values_missing_values(self):
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
    def test_divide_values_weird_value(self):
        """tests if a non float or integer value is encountered
        should raise an error"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered")
        fp = open(fpath, "w")
        fp.write("	sample1	sample2\n")
        fp.write("Cluster0	30.2	15.2\n")
        fp.write("Cluster1	40.5	ABCDE")
        fp.close()
        self.assertRaises(TypeError, divide_values, {'Cluster2': '60.6', 'Cluster0': '30.2', 'Cluster1': '40.5'})
        shutil.rmtree(tdir)

class Test9(unittest.TestCase):
    def test_translate_genes_basic_function(self):
        """tests to see if the translation is correct, and if shorter sequences
        get filtered out"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered")
        fp = open(fpath, "w")
        fp.write(">gi|22123922|ref|NC_004088.1|_3285\n")
        fp.write("ATGAATCCTCACCTAACCGAACACCCCCCAGTCGGGGATATTGACGCCCTGTTGCAGGACACCTGGCTACAGGTGATCAGCCTGCGTCAAGGGGTAACCTGTGCCGAGGGCGAAGGGCAGGCATTCTGGCAGCGCTGTGTGGCGGACATTGAACGTGTCCATCAGGCGCTGAAAGACGCCGGTCACAGCGAGCAGAGTTGCCAGCACATCCGATACGCCCAATGTGCACTGCTGGATGAG\n")
        fp.write(">gi|22123922|ref|NC_004088.1|_1575\n")
        fp.write("ATGAAGCTAAATATCAAAGTTAATTGTTCTTATATCTGTGAACCCATACGTAAGCAA")
        fp.close()
        self.assertEqual(translate_genes(fpath, "tmp.out", 10), 'MNPHLTEHPPVGDIDALLQDTWLQVISLRQGVTCAEGEGQAFWQRCVADIERVHQALKDAGHSEQSCQHIRYAQCALLDE')
        shutil.rmtree(tdir)
    def test_translate_genes_out_of_frame(self):
        """test the condition where the sequence is not in frame"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered")
        fp = open(fpath, "w")
        fp.write(">gi|22123922|ref|NC_004088.1|_3285\n")
        fp.write("CTCATCCAGCAGTGCACATTGGGCGTATCGGATGTGCTGGCAACTCTGCTCGCTGTGACCGGCGTCTTTCAGCGCCTGATGGACACGTTCAATGTCCGCCACACAGCGCTGCCAGAATGCCTGCCCTTCGCCCTCGGCACAGGTTACCCCTTGACGCAGGCTGATCACCTGTAGCCAGGTGTCCTGCAACAGGGCGTCAATATCCCCGACTGGGGGGTGTTCGGTTAGGTGAGGATTCAT")
        fp.close()
        self.assertEqual(translate_genes(fpath, "tmp.out", 10), 'LIQQCTLGVSDVLATLLAVTGVFQRLMDTFNVRHTALPECLPFALGTGYPLTQADHL')
        shutil.rmtree(tdir)
    def test_translate_genes_odd_characters(self):
        """tests the condition where weird characters are encountered"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered")
        fp = open(fpath, "w")
        fp.write(">gi|22123922|ref|NC_004088.1|_3285\n")
        fp.write("1234567890")
        fp.close()
        self.assertRaises(TypeError, translate_genes, fpath)
        shutil.rmtree(tdir)
    def test_translate_genes_non_fasta(self):
        """tests the condition where the file is not fasta.  Will
        not throw an error, but will report an empty set"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered")
        fp = open(fpath, "w")
        fp.write("not a fasta file")
        fp.close()
        self.assertEqual(translate_genes(fpath, "tmp.out", 10), None)
        shutil.rmtree(tdir)

class Test10(unittest.TestCase):
    def test_rename_fasta_basic_function(self):
        """tests standard functionality of rename_fasta_header function"""
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
    def test_rename_fasta_odd_characters(self):
        """tests condition where non-normal characters are encountered"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered")
        fp = open(fpath, "w")
        fp.write(">gi|22123922|ref|NC_004088.1|_3285\n")
        fp.write("1234567890")
        fp.close()
        self.assertEqual(rename_fasta_header(fpath, "tmp.out"), ['>centroid_gi|22123922|ref|NC_004088.1|_3285'])
        os.system("rm tmp.out")
        shutil.rmtree(tdir)
    def test_rename_fasta_non_fasta(self):
        """tests condition where the file is not fasta"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered")
        fp = open(fpath, "w")
        fp.write("not a fasta file")
        fp.close()
        self.assertEqual(rename_fasta_header(fpath, "tmp.out"), [])
        os.system("rm tmp.out")
        shutil.rmtree(tdir)

class Test11(unittest.TestCase):
    def test_auto_increment_function(self):
        self.assertEqual(autoIncrement(), 5)

class Test12(unittest.TestCase):
    def test_prune_matrix_basic_function(self):
        """tests the basic functionaly of the prune_matrix function"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered")
        fp = open(fpath, "w")
        fp.write("        E2348_69_all    H10407_all      O157_H7_sakai_all       SSON_046_all\n")
        fp.write("IpaH3   0.03    0.03    0.03    1.00\n")
        fp.write("LT      0.00    1.00    0.00    0.00\n")
        fp.write("ST1     0.00    1.00    0.12    0.12\n")
        fp.write("bfpB    1.00    0.00    0.00    0.00\n")
        fp.write("stx2a   0.07    0.08    0.98    0.07\n")
        fp.close()
        npath = os.path.join(tdir,"group1")
        np = open(npath, "w")
        np.write("E2348_69_all")
        np.close()
        opath = os.path.join(tdir,"group2")
        op = open(opath, "w")
        op.write("H10407_all")
        op.close()
        self.assertEqual(prune_matrix(fpath,npath,opath), (['E2348_69_all'], ['H10407_all'], [0, 2, 3, 4], [0, 1, 3, 4]))
        shutil.rmtree(tdir)
        os.system("rm group*_pruned.txt")
    def test_prune_matrix_empty_file(self):
        """tests the case where one of the group files is empty"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered")
        fp = open(fpath, "w")
        fp.write("        E2348_69_all    H10407_all      O157_H7_sakai_all       SSON_046_all\n")
        fp.write("IpaH3   0.03    0.03    0.03    1.00\n")
        fp.write("LT      0.00    1.00    0.00    0.00\n")
        fp.write("ST1     0.00    1.00    0.12    0.12\n")
        fp.write("bfpB    1.00    0.00    0.00    0.00\n")
        fp.write("stx2a   0.07    0.08    0.98    0.07\n")
        fp.close()
        npath = os.path.join(tdir,"group1")
        np = open(npath, "w")
        np.write("E2348_69_all")
        np.close()
        opath = os.path.join(tdir,"group2")
        op = open(opath, "w")
        op.write("")
        op.close()
        self.assertEqual(prune_matrix(fpath,npath,opath), (['E2348_69_all'], [], [0, 2, 3, 4], [0, 1, 2, 3, 4]))
        shutil.rmtree(tdir)
        os.system("rm group*_pruned.txt")
    def test_prune_matrix_no_match(self):
        """tests the case where the genome file isn't in the matrix"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered")
        fp = open(fpath, "w")
        fp.write("        E2348_69_all    H10407_all      O157_H7_sakai_all       SSON_046_all\n")
        fp.write("IpaH3   0.03    0.03    0.03    1.00\n")
        fp.write("LT      0.00    1.00    0.00    0.00\n")
        fp.write("ST1     0.00    1.00    0.12    0.12\n")
        fp.write("bfpB    1.00    0.00    0.00    0.00\n")
        fp.write("stx2a   0.07    0.08    0.98    0.07\n")
        fp.close()
        npath = os.path.join(tdir,"group1")
        np = open(npath, "w")
        np.write("not_there")
        np.close()
        opath = os.path.join(tdir,"group2")
        op = open(opath, "w")
        op.write("not this one either")
        op.close()
        self.assertEqual(prune_matrix(fpath,npath,opath), (['not_there'], ['not this one either'], [0, 1, 2, 3, 4], [0, 1, 2, 3, 4]))
        shutil.rmtree(tdir)
        os.system("rm group*_pruned.txt")

class Test13(unittest.TestCase):
    def test_compare_values_basic_function(self):
        """test basic functionality of compare_values function.  Tests
        that the values are being parsed, and that the Numpy mean feature
        is working correctly"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"group1_pruned")
        fp = open(fpath, "w")
        fp.write("        E2348_69_all    H10407_all\n")
        fp.write("IpaH3   0.03    0.03\n")
        fp.write("LT      0.00    1.00\n")
        fp.write("ST2     0.00    0.92\n")
        fp.write("bfpB    1.00    0.00\n")
        fp.write("stx2a   0.07    0.08")
        fp.close()
        npath = os.path.join(tdir,"group2_pruned")
        np = open(npath, "w")
        np.write("        SSON_046_all\n")
        np.write("IpaH3   0.03\n")
        np.write("LT      1.00\n")
        np.write("ST2     1.00\n")
        np.write("bfpB    0.00\n")
        np.write("stx2a   0.08")
        np.close()
        self.assertEqual(compare_values(fpath, npath, "0.8", "0.4"), ([1.00, 0.92, 1.00], [1.00, 1.00], [0.03, 0.5, 0.46, 0.5, 0.07500000000000001]))
        shutil.rmtree(tdir)
    def test_compare_values_border_cases(self):
        """tests the condition where BSR values are near the border regions
        differentiated by the function"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"group1_pruned")
        fp = open(fpath, "w")
        fp.write("		E2348_69_all\n")
        fp.write("IpaH3 	0.03\n")
        fp.write("LT 	0.00\n")
        fp.write("ST2 	0.00\n")
        fp.write("bfpB 	0.81\n")
        fp.write("stx2a 	0.07")
        fp.close()
        npath = os.path.join(tdir,"group2_pruned")
        np = open(npath, "w")
        np.write("        H10407_all\n")
        np.write("IpaH3   0.03\n")
        np.write("LT      0.80\n")
        np.write("ST2     1.00\n")
        np.write("bfpB    0.00\n")
        np.write("stx2a   0.79")
        np.close()
        self.assertEqual(compare_values(fpath, npath, "0.8", "0.4"), ([0.81], [0.80, 1.00], [0.03, 0.0, 0.0, 0.81, 0.07]))
        shutil.rmtree(tdir)
        os.system("rm group*_out.txt")

class Test14(unittest.TestCase):
    def test_find_uniques_basic_function(self):
        """tests the basic functionality of the find_uniques function"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"combined.txt")
        fp = open(fpath, "w")
        fp.write("IpaH3   0.03    0       1       0       0.03    0       1       0\n")
        fp.write("LT 	0.0 	0 	1 	0	0.8 	1 	1 	1\n")
        fp.write("ST2 	0.0 	0 	1 	0	1.0 	1 	1 	1\n")
        fp.write("bfpB 	0.81 	1 	1 	1	0.0 	0 	1 	0\n")
        fp.write("stx2a 	0.07 	0 	1 	0	0.79 	0 	1 	1")
        fp.close()
        npath = os.path.join(tdir,"fasta")
        np = open(npath, "w")
        np.write(">bfpB\n")
        np.write("ATGAAACTTGGCAGGTATTCACTTTTCTTATTG\n")
        np.write(">LT\n")
        np.write("ATGCCCAGAGGGCATAATGAGTACTTCGA\n")
        np.write(">ST2\n")
        np.write("ATGAAGAAATCAATATTATTTATTTTTCTTTCTGTATTGTCTTTT")
        np.close()
        self.assertEqual(find_uniques(fpath, npath), (['bfpB'], ['LT', 'ST2'], ['bfpB']))
        shutil.rmtree(tdir)
        os.system("rm group*_unique_seqs.fasta")
    def test_find_uniques_mismatched_headers(self):
        """tests the case where headers don't match to those headers
        in a matrix"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"combined.txt")
        fp = open(fpath, "w")
        fp.write("IpaH3\t0.03\t0\t1\t0\t0.03\t0\t1\t0\n")
        fp.write("LT 	0.0 	0 	1 	0	0.8 	1 	1 	1\n")
        fp.write("ST2 	0.0 	0 	1 	0	1.0 	1 	1 	1\n")
        fp.write("bfpB 	0.81 	1 	1 	1	0.0 	0 	1 	0\n")
        fp.write("stx2a 	0.07 	0 	1 	0	0.79 	0 	1 	1")
        fp.close()
        npath = os.path.join(tdir,"fasta")
        np = open(npath, "w")
        np.write(">nomatch\n")
        np.write("ATGAAACTTGGCAGGTATTCACTTTTCTTATTG\n")
        np.write(">nomatch1\n")
        np.write("ATGCCCAGAGGGCATAATGAGTACTTCGA\n")
        np.write(">nomatch2\n")
        np.write("ATGAAGAAATCAATATTATTTATTTTTCTTTCTGTATTGTCTTTT")
        np.close()
        self.assertEqual(find_uniques(fpath, npath), (['bfpB'], ['LT', 'ST2'], []))
        shutil.rmtree(tdir)
        os.system("rm group*_unique_seqs.fasta")

class Test15(unittest.TestCase):
    def test_filter_genomes_basic_function(self):
        """test the basic functionality of the filter_genomes function"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("        E2348_69_all    H10407_all      O157_H7_sakai_all       SSON_046_all\n")
        fp.write("IpaH3   0.03    0.03    0.03    1.00\n")
        fp.write("LT      0.00    1.00    0.00    0.00\n")
        fp.write("ST1     0.00    1.00    0.12    0.12\n")
        fp.write("bfpB    1.00    0.00    0.00    0.00\n")
        fp.write("stx2a   0.07    0.08    0.98    0.07\n")
        fp.close()
        npath = os.path.join(tdir,"genomes")
        np = open(npath, "w")
        np.write("H10407_all\n")
        np.write("SSON_046_all")
        np.close()
        """1 and 3 should be the correct indices, based on the input data.  These
        are the genomes that should be removed"""
        self.assertEqual(filter_genomes(npath, fpath), [1, 3])
        shutil.rmtree(tdir)
    def test_filter_genomes_no_matches(self):
        """tests the condition where the group names provided are not in
        the input matrix"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("        E2348_69_all    H10407_all      O157_H7_sakai_all       SSON_046_all\n")
        fp.write("IpaH3   0.03    0.03    0.03    1.00\n")
        fp.write("LT      0.00    1.00    0.00    0.00\n")
        fp.write("ST1     0.00    1.00    0.12    0.12\n")
        fp.write("bfpB    1.00    0.00    0.00    0.00\n")
        fp.write("stx2a   0.07    0.08    0.98    0.07\n")
        fp.close()
        npath = os.path.join(tdir,"genomes")
        np = open(npath, "w")
        np.write("absent\n")
        np.write("absent2")
        np.close()
        self.assertEqual(filter_genomes(npath, fpath), [])
        shutil.rmtree(tdir)

class Test16(unittest.TestCase):
    def test_filter_matrix_basic_function(self):
        """test the basic functionality of the filter_matrix function"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("        E2348_69_all    H10407_all      O157_H7_sakai_all       SSON_046_all\n")
        fp.write("IpaH3   0.03    0.03    0.03    1.00\n")
        fp.write("LT      0.00    1.00    0.00    0.00\n")
        fp.write("ST1     0.00    1.00    0.12    0.12\n")
        fp.write("bfpB    1.00    0.00    0.00    0.00\n")
        fp.write("stx2a   0.07    0.08    0.98    0.07\n")
        fp.close()
        self.assertEqual(filter_matrix([1, 3], fpath, "test"), [['', 'E2348_69_all', 'O157_H7_sakai_all'],['IpaH3', '0.03', '0.03'],['LT', '0.00', '0.00'],['ST1', '0.00', '0.12'],['bfpB', '1.00', '0.00'],['stx2a', '0.07', '0.98']])
        shutil.rmtree(tdir)
        os.system("rm test_genomes.matrix")
    def test_filter_matrix_non_numerics(self):
        """tests the case where non-numeric values are in matrix"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("        E2348_69_all    H10407_all      O157_H7_sakai_all       SSON_046_all\n")
        fp.write("IpaH3   cats    0.03    dogs    1.00\n")
        fp.write("LT      0.00    1.00    0.00    0.00\n")
        fp.write("ST1     0.00    1.00    0.12    0.12\n")
        fp.write("bfpB    1.00    0.00    0.00    0.00\n")
        fp.write("stx2a   0.07    0.08    0.98    0.07\n")
        fp.close()
        self.assertEqual(filter_matrix([1, 3], fpath, "test"), [['', 'E2348_69_all', 'O157_H7_sakai_all'],['IpaH3', 'cats', 'dogs'],['LT', '0.00', '0.00'],['ST1', '0.00', '0.12'],['bfpB', '1.00', '0.00'],['stx2a', '0.07', '0.98']])
        shutil.rmtree(tdir)
        os.system("rm test_genomes.matrix")

class Test17(unittest.TestCase):
    def test_get_core_gene_stats_basic_function(self):
        """test the basic functionality of the get_core_gene_stats function"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("        E2348_69_all    H10407_all      O157_H7_sakai_all       SSON_046_all\n")
        fp.write("IpaH3   1.00    1.00    1.00    1.00\n")
        fp.write("LT      0.00    1.00    0.00    0.00\n")
        fp.write("ST1     0.00    1.00    0.12    0.12\n")
        fp.write("bfpB    1.00    0.00    0.00    0.00\n")
        fp.write("stx2a   0.07    0.08    0.98    0.07\n")
        fp.close()
        self.assertEqual(get_core_gene_stats(fpath, 0.8, 0.4, 0), (1, 4))
        shutil.rmtree(tdir)
        os.system("rm core_gene_ids.txt unique_gene_ids.txt")
    def test_get_core_gene_stats_empty_line(self):
        """tests the case where an empty line is found"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("        E2348_69_all    H10407_all      O157_H7_sakai_all       SSON_046_all\n")
        fp.write("IpaH3   1.00    1.00    1.00    1.00\n")
        fp.write("LT      0.00    1.00    0.00    0.00\n")
        fp.write("ST1     0.00    1.00    0.12    0.12\n")
        fp.write("bfpB    1.00    0.00    0.00    0.00\n")
        fp.write("stx2a\n")
        fp.close()
        self.assertRaises(TypeError, get_core_gene_stats, fpath, 0.8, 0.4, 0)
        shutil.rmtree(tdir)
        os.system("rm core_gene_ids.txt unique_gene_ids.txt")

class Test18(unittest.TestCase):
    def test_get_frequencies_basic_function(self):
        """test the basic functionality of the get_frequencies function"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("        E2348_69_all    H10407_all      O157_H7_sakai_all       SSON_046_all\n")
        fp.write("IpaH3   1.00    1.00    1.00    1.00\n")
        fp.write("LT      0.00    1.00    0.00    0.00\n")
        fp.write("ST1     0.00    1.00    0.12    0.12\n")
        fp.write("bfpB    1.00    0.00    0.00    0.00\n")
        fp.write("stx2a   0.07    0.08    0.98    0.07\n")
        fp.close()
        self.assertEqual(get_frequencies(fpath, 0.8), {1:4,4:1})
        shutil.rmtree(tdir)
    def test_get_frequencies_borders(self):
        """test the border case"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("        E2348_69_all    H10407_all      O157_H7_sakai_all       SSON_046_all\n")
        fp.write("IpaH3   0.80    1.00    1.00    1.00\n")
        fp.write("LT      0.00    1.00    0.79    0.00\n")
        fp.write("ST1     0.00    1.00    0.12    0.12\n")
        fp.write("bfpB    1.00    0.00    0.00    0.00\n")
        fp.write("stx2a   0.07    0.08    0.98    0.07\n")
        fp.close()
        self.assertEqual(get_frequencies(fpath, 0.8), {1:4,4:1})
        shutil.rmtree(tdir)
        os.system("rm frequency_data.txt")

class Test20(unittest.TestCase):
    def test_filter_paralogs_basic_function(self):
        """tests the basic functionality of the filter_paralogs function"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("        E2348_69_all    H10407_all      O157_H7_sakai_all       SSON_046_all\n")
        fp.write("IpaH3   0.80    1.00    1.00    1.00\n")
        fp.write("LT      0.00    1.00    0.79    0.00\n")
        fp.write("ST1     0.00    1.00    0.12    0.12\n")
        fp.write("bfpB    1.00    0.00    0.00    0.00\n")
        fp.write("stx2a   0.07    0.08    0.98    0.07\n")
        fp.close()
        npath = os.path.join(tdir,"name_file")
        np = open(npath, "w")
        np.write("LT\n")
        np.write("stx2a")
        np.close()
        #self.assertEqual(filter_paralogs(fpath, npath), ['IpaH3', 'ST1', 'bfpB'])
        self.assertEqual(filter_paralogs(fpath, npath), 2)
        shutil.rmtree(tdir)
        os.system("rm bsr_matrix_values_filtered.txt")
    def test_filter_paralogs_no_matches(self):
        """tests the case where no matching names are found"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("        E2348_69_all    H10407_all      O157_H7_sakai_all       SSON_046_all\n")
        fp.write("IpaH3   0.80    1.00    1.00    1.00\n")
        fp.write("LT      0.00    1.00    0.79    0.00\n")
        fp.write("ST1     0.00    1.00    0.12    0.12\n")
        fp.write("bfpB    1.00    0.00    0.00    0.00\n")
        fp.write("stx2a   0.07    0.08    0.98    0.07\n")
        fp.close()
        npath = os.path.join(tdir,"name_file")
        np = open(npath, "w")
        np.write("no_match")
        np.close()
        #self.assertEqual(filter_paralogs(fpath, npath), ['IpaH3', 'LT', 'ST1', 'bfpB', 'stx2a'])
        self.assertEqual(filter_paralogs(fpath, npath), 0)
        shutil.rmtree(tdir)
        os.system("rm bsr_matrix_values_filtered.txt")

class Test21(unittest.TestCase):
    def test_filter_variome_basic_function(self):
        """tests the basic functionality of the filter_variome function"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("        E2348_69_all    H10407_all      O157_H7_sakai_all       SSON_046_all\n")
        fp.write("IpaH3   0.80    1.00    1.00    1.00\n")
        fp.write("LT      0.00    1.00    0.79    0.00\n")
        fp.write("ST1     0.00    1.00    0.12    0.12\n")
        fp.write("bfpB    1.00    0.00    0.00    0.00\n")
        fp.write("stx2a   0.07    0.08    0.98    0.07\n")
        fp.close()
        self.assertEqual(filter_variome(fpath, "0.8", "1"), ['LT', 'ST1', 'bfpB', 'stx2a'])
        shutil.rmtree(tdir)
        os.system("rm variome_BSR_matrix")
    def test_filter_variome_bad_input(self):
        """tests the case where a malformed file is found"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("        E2348_69_all    H10407_all      O157_H7_sakai_all       SSON_046_all\n")
        fp.write("IpaH3   0.80    1.00    1.00    1.00\n")
        fp.write("LT      0.00    1.00    0.79    0.00\n")
        fp.write("ST1     0.00    1.00    0.12    0.12\n")
        fp.write("bfpB    1.00    0.00    0.00    0.00\n")
        fp.write("")
        fp.close()
        self.assertEqual(filter_variome(fpath, "0.8", "1"), ['LT', 'ST1', 'bfpB'])
        shutil.rmtree(tdir)
        os.system("rm variome_BSR_matrix")

class Test22(unittest.TestCase):
    def test_process_pangenome_basic_function(self):
        """tests basic functionality of the process_pangenome function"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        print(tdir)
        fp = open(fpath, "w")
        fp.write("\tSSON_046_all\n")
        fp.write("IpaH3\t1.00\n")
        fp.write("LT\t0.00\n")
        fp.write("ST1\t0.12\n")
        fp.write("bfpB\t0.00\n")
        fp.write("stx2a\t0.07\n")
        fp.write("junk\t1.00\n")
        fp.close()
        self.assertEqual(process_pangenome(fpath, "0.8", "0.4", 1, "all","random"), ([2.0], [2.0], [2.0]))
        #shutil.rmtree(tdir)
        #os.system("rm core_replicates.txt uniques_replicates.txt accumulation_replicates.txt")
    def test_process_pangenome_cores_only(self):
        """tests when you only run the cores"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("\tSSON_046_all\n")
        fp.write("IpaH3\t1.00\n")
        fp.write("LT\t0.00\n")
        fp.write("ST1\t0.12\n")
        fp.write("bfpB\t0.00\n")
        fp.write("stx2a\t0.07\n")
        fp.write("junk\t1.00\n")
        fp.close()
        self.assertEqual(process_pangenome(fpath, "0.8", "0.4", 1, "core","random"), ([],[],[2.0]))
        shutil.rmtree(tdir)
        #os.system("rm core_replicates.txt")
    def test_process_pangenome_unis_only(self):
        """tests when you only want the uniques"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("\tSSON_046_all\n")
        fp.write("IpaH3\t1.00\n")
        fp.write("LT\t0.00\n")
        fp.write("ST1\t0.12\n")
        fp.write("bfpB\t0.00\n")
        fp.write("stx2a\t0.07\n")
        fp.write("junk\t1.00\n")
        fp.close()
        self.assertEqual(process_pangenome(fpath, "0.8", "0.4", 1, "uni","random"), ([],[2.0],[]))
        shutil.rmtree(tdir)
        #os.system("rm uniques_replicates.txt")
    def test_process_pangenome_core_only(self):
        """tests when you only want the accums"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("\tSSON_046_all\n")
        fp.write("IpaH3\t1.00\n")
        fp.write("LT\t0.00\n")
        fp.write("ST1\t0.12\n")
        fp.write("bfpB\t0.00\n")
        fp.write("stx2a\t0.07\n")
        fp.write("junk\t1.00\n")
        fp.close()
        self.assertEqual(process_pangenome(fpath, "0.8", "0.4", 1, "acc", "random"), ([2.0],[],[]))
        shutil.rmtree(tdir)
        #os.system("rm accumulation_replicates.txt")
    def test_process_pangenome_upper(self):
        """tests when you only want the accums"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("\tSSON_046_all\n")
        fp.write("IpaH3\t0.89\n")
        fp.write("LT\t0.00\n")
        fp.close()
        self.assertEqual(process_pangenome(fpath, "0.9", "0.4", 1, "all", "random"), ([1.0], [1.0], [1.0]))
        shutil.rmtree(tdir)
        #os.system("rm core_replicates.txt uniques_replicates.txt accumulation_replicates.txt")
    def test_process_pangenome_upper(self):
        """tests when you only want the accums"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("\tSSON_046_all\n")
        fp.write("IpaH3\t0.91\n")
        fp.write("LT\t0.41\n")
        fp.close()
        self.assertEqual(process_pangenome(fpath, "0.9", "0.4", 1, "all", "random"), ([1.0], [1.0], [1.0]))
        shutil.rmtree(tdir)
        #os.system("rm random_core_replicates.txt random_uniques_replicates.txt random_accumulation_replicates.txt")

class Test23(unittest.TestCase):
    def test_bsr_to_pangb_basic_function(self):
        """tests the basic functionality"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("        E2348_69_all    H10407_all      O157_H7_sakai_all       SSON_046_all\n")
        fp.write("ST1     0.00    1.00    0.12    0.12\n")
        fp.write("bfpB    1.00    0.81    0.00    0.00\n")
        fp.write("")
        fp.close()
        self.assertEqual(bsr_to_pangp(fpath,0.8), ['bfpB','1','1','-','-'])
        shutil.rmtree(tdir)
        os.system("rm panGP_matrix.txt")
    def test_bsr_to_pangb_equal_function(self):
        """tests the basic functionality"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"sample_matrix.txt")
        fp = open(fpath, "w")
        fp.write("        E2348_69_all    H10407_all      O157_H7_sakai_all       SSON_046_all\n")
        fp.write("ST1     0.00    1.00    0.12    0.12\n")
        fp.write("bfpB    1.00    0.80    0.00    0.00\n")
        fp.write("")
        fp.close()
        self.assertEqual(bsr_to_pangp(fpath,0.8), ['bfpB','1','1','-','-'])
        shutil.rmtree(tdir)
        os.system("rm panGP_matrix.txt")

class Test24(unittest.TestCase):
    def test_get_cluster_ids_basic_function(self):
        """basic functionality"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        npath = os.path.join(tdir,"fasta")
        np = open(npath, "w")
        np.write(">bfpB\n")
        np.write("ATGAAACTTGGCAGGTATTCACTTTTCTTATTG\n")
        np.write(">LT\n")
        np.write("ATGCCCAGAGGGCATAATGAGTACTTCGA\n")
        np.write(">ST2\n")
        np.write("ATGAAGAAATCAATATTATTTATTTTTCTTTCTGTATTGTCTTTT")
        np.close()
        self.assertEqual(get_cluster_ids(npath), ['bfpB','LT','ST2'])
        shutil.rmtree(tdir)
    def test_get_cluster_ids_weird_characters(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        npath = os.path.join(tdir,"fasta")
        np = open(npath, "w")
        np.write(">bfp-B\n")
        np.write("ATGAAACTTGGCAGGTATTCACTTTTCTTATTG\n")
        np.write(">LT_X\n")
        np.write("ATGCCCAGAGGGCATAATGAGTACTTCGA\n")
        np.write(">ST#%$\n")
        np.write("ATGAAGAAATCAATATTATTTATTTTTCTTTCTGTATTGTCTTTT")
        np.close()
        self.assertEqual(get_cluster_ids(npath), ['bfp-B','LT_X','ST#%$'])

class Test25(unittest.TestCase):
    def test_intergenics_basic_function(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        npath = os.path.join(tdir,"fasta")
        np = open(npath, "w")
        np.write(">test\n")
        np.write("ATGAAACTTGGCAGGTATTCACTTTTCTTATTGATGAAACTTGGCAGGTATTCACTTTTCTTATTGATGAAACTTGGCAGGTATTCACTTTTCTTATTG\n")
        np.close()
        npath2 = os.path.join(tdir,"ranges")
        np = open(npath2, "w")
        np.write("test\t2\t54\n")
        np.close()
        self.assertEqual(parse_ranges_file(npath,npath2,"test","true"), ['TGAAACTTGGCAGGTATTCACTTTTCTTATTGATGAAACTTGGCAGGTATTCA'])
        os.system("rm test.intergenics.seqs")
        shutil.rmtree(tdir)

if __name__ == "__main__":
    unittest.main()
    main()
