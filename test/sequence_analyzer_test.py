import unittest
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from codon_tools.sequence_analyzer import SequenceAnalyzer
from codon_tools.lookup_tables import opt_codons_E_coli, reverse_genetic_code

class TestSequenceAnalyzer(unittest.TestCase):
    def setUp(self):
        self.sa = SequenceAnalyzer()
        # seq1 is wt GFP sequence
        self.seq1 = Seq('ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG', IUPAC.unambiguous_dna)
        self.seq1_codon_counts = {'GCA': 0, 'GCC': 8, 'GCG': 0, 'GCT': 0,
            'AGA': 0, 'AGG': 0, 'CGA': 0, 'CGT': 0, 'CGC': 6, 'CGG': 0, 
            'AAC': 13, 'AAT': 0, 'GAC': 16, 'GAT': 2, 'TGC': 2, 'TGT': 0,
            'CAA': 0, 'CAG': 8, 'GAA': 1, 'GAG': 15, 'GGA': 0, 'GGC': 19,
            'GGG': 3, 'GGT': 0, 'CAC': 9, 'CAT': 0, 'ATA': 0, 'ATC': 12, 
            'ATT': 0, 'CTA': 0, 'CTC': 3, 'CTG': 18, 'CTT': 0, 'TTA': 0,
            'TTG': 0, 'TTT': 0, 'TTC': 12, 'CCA': 0, 'CCC': 10, 'CCG': 0,
            'CCT': 0, 'AGC': 7, 'AGT': 0, 'TCA': 0, 'TCC': 3, 'TCG': 0,
            'TCT': 0, 'ACA': 0, 'ACT': 1, 'ACC': 15, 'ACG': 0, 'TAC': 10,
            'TAT': 1, 'GTA': 1, 'GTC': 4, 'GTG': 13, 'GTT': 0, 'TGG': 1,
            'ATG': 6, 'AAA': 1, 'AAG': 19}
    
    def test_count_muts_to_stop(self):
        m2stop_count = self.sa.count_muts_to_stop(self.seq1)
        self.assertEqual(m2stop_count, 70)
    
    def test_count_CpG(self):
        CpG_count, UpA_count = self.sa.count_CpG(self.seq1)
        self.assertEqual(CpG_count, 60)
        self.assertEqual(UpA_count, 13)

    def test_calc_Fop(self):
        opt_count, total_count, Fop = self.sa.calc_Fop(self.seq1, opt_codons = opt_codons_E_coli)
        self.assertEqual(opt_count, 146)
        self.assertEqual(total_count, 212)
        self.assertEqual(Fop, opt_count/total_count)
        
    def test_count_diffs(self):
        s1 = 'AGCT'
        s2 = 'GGCT'
        s3 = 'AGCA'
        s4 = 'GC-A'
        self.assertEqual(self.sa.count_diffs(s1, s2), 1)
        self.assertEqual(self.sa.count_diffs(s1, s3), 1)
        self.assertEqual(self.sa.count_diffs(s2, s3), 2)
        self.assertEqual(self.sa.count_diffs(s3, s3), 0)
        self.assertEqual(self.sa.count_diffs(s1, s4), 4)
        self.assertEqual(self.sa.count_diffs(s2, s4), 3)
        
    def test_count_codons(self):
        counts = self.sa.count_codons(self.seq1)
        for codon, count in counts.items():
            self.assertEqual(count, self.seq1_codon_counts[codon])

    def test_calc_syn_codon_freqs(self):
        freqs = self.sa.calc_syn_codon_freqs(self.seq1)

        # calculate frequencies from stored codon counts 
        counts = self.seq1_codon_counts.copy()       
        for aa in reverse_genetic_code:
            if aa == '*': # we skip over stop codons
                continue
            codons = reverse_genetic_code[aa]
            fam_counts = {c:counts[c] for c in codons}
            total = sum(fam_counts.values())
            if total > 0:
                for c in fam_counts:
                    counts[c] = fam_counts[c]/total
        # now compare
        for codon, freq in freqs.items():
            self.assertEqual(freq, counts[codon])
    
    def tearDown(self):
        # Do any clean up here
        pass
            


