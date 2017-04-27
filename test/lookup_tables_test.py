import unittest
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from codon_tools.lookup_tables import opt_codons_E_coli, reverse_genetic_code

class TestLookupTables(unittest.TestCase):
    def test_reverse_genetic_code(self):
        tested_codons = {}
        for aa, codons in reverse_genetic_code.items():
            for codon in codons:
                self.assertEqual(aa, Seq(codon, generic_dna).translate())
                if codon in tested_codons:
                    self.assertTrue(False)
                else:
                    tested_codons[codon] = 1

    def test_opt_codons_E_coli(self):
        # known optimal codons for E. coli
        opt_codons = { 'A':['GCT'], 'R':['CGT', 'CGC'], 'N':['AAC'], 'D':['GAC'], 'C':['TGC'], 'Q':['CAG'], 'E':['GAA'], 'G':['GGT','GGC'], 'H':['CAC'], 'I':['ATC'], 'L':['CTG'], 'F':['TTC'], 'P':['CCG'], 'S':['TCT','TCC'], 'T':['ACT','ACC'], 'Y':['TAC'], 'V':['GTT','GTA'] }
        
        for aa, codons in opt_codons_E_coli.items():
            observed = set(codons)
            expected = set(opt_codons[aa])
            self.assertEqual(observed, expected)
            

