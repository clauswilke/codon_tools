import unittest
from codon_tools.codon_analyzer import CodonAnalyzer

class TestCodonAnalyzer(unittest.TestCase):
    def setUp(self):
        # Do any set up that you need here
        self.my_analyzer = CodonAnalyzer()
        self.mut_to_stop = { 'AAA': 1, 'AAG': 1, 'AAC': 0, 'AAT': 0,
         'AGA': 1, 'AGG': 0, 'AGC': 0, 'AGT': 0, 'ACA': 0, 'ACG': 0,
         'ACC': 0, 'ACT': 0, 'ATA': 0, 'ATG': 0, 'ATC': 0, 'ATT': 0,
         'GAA': 1, 'GAG': 1, 'GAC': 0, 'GAT': 0, 'GGA': 1, 'GGG': 0,
         'GGC': 0, 'GGT': 0, 'GCA': 0, 'GCG': 0, 'GCC': 0, 'GCT': 0,
         'GTA': 0, 'GTG': 0, 'GTC': 0, 'GTT': 0, 'CAA': 1, 'CAG': 1,
         'CAC': 0, 'CAT': 0, 'CGA': 1, 'CGG': 0, 'CGC': 0, 'CGT': 0,
         'CCA': 0, 'CCG': 0, 'CCC': 0, 'CCT': 0, 'CTA': 0, 'CTG': 0,
         'CTC': 0, 'CTT': 0, 'TAA': 2, 'TAG': 1, 'TAC': 2, 'TAT': 2,
         'TGA': 1, 'TGG': 2, 'TGC': 1, 'TGT': 1, 'TCA': 2, 'TCG': 1,
         'TCC': 0, 'TCT': 0, 'TTA': 2, 'TTG': 1, 'TTC': 0, 'TTT': 0 }

    def test_count_mutations_to_stop(self):
        # loop over all codons and compare calculated result to stored result
        for seq, count in self.mut_to_stop.items():
            self.assertEqual(self.my_analyzer.count_stop_muts(seq), count)
