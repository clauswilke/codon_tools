import random

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

reverse_genetic_code = {  'A':['GCA', 'GCC', 'GCG', 'GCT'], 'R':['AGA', 'AGG', 'CGA', 'CGT', 'CGC', 'CGG'], 'N':['AAC', 'AAT'], 'D':['GAC', 'GAT'], 'C':['TGC', 'TGT' ], 'Q':['CAA', 'CAG'], 'E':['GAA', 'GAG'], 'G':['GGA', 'GGC', 'GGG', 'GGT'], 'H':['CAC', 'CAT'], 'I':['ATA', 'ATC', 'ATT'], 'L':['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'], 'F':['TTT', 'TTC'], 'P':['CCA', 'CCC', 'CCG', 'CCT'], 'S':['AGC', 'AGT', 'TCA', 'TCC','TCG','TCT'], 'T':[ 'ACA', 'ACT', 'ACC', 'ACG' ], 'Y':['TAC', 'TAT'], 'V':['GTA', 'GTC', 'GTG', 'GTT'], 'W':['TGG'], 'M':['ATG'], 'K':['AAA', 'AAG']}

class CodonOptimizer:
    """Class to optimize codon usage in an ORF. When instantiated, it needs to be given a ``scorer`` object that implements a ``score`` function. 
    """
    def __init__(self, scorer):
        self.scorer = scorer

    def random_reverse_translate(self, protein):
        # choose random codon for each amino acid
        dna = ''.join([random.choice(reverse_genetic_code[aa]) for aa in protein.upper() ])
        return Seq(dna, IUPAC.unambiguous_dna)

    def change_random_codon(self, seq, start_window, end_window):
        """``seq``: The nucleotide sequence in which a codon should be changed
``start_window``: number of codons to omit at beginning of sequence
``end_window``: number of codons to omit at end of sequence
"""
        #print len(seq), str(seq)
        assert len( seq ) % 3 == 0
    
        start = start_window
        end = len( seq )/3 - end_window - 1
        assert start <= end
        aa_index = random.randint(start, end)

        nt_index = 3*aa_index
    
        return seq[0:nt_index] + \
          self.random_reverse_translate(seq[nt_index:nt_index+3].translate()) + \
          seq[nt_index+3:]

    def hillclimb(self, seq, start_window = 0, end_window = 0, target_score = None, tolerance = .01, max_wait_count = 1000, maximize = True, verbosity = 1 ):
        '''The hillclimb will stop if the obtained score is within ``tolerance`` of ``target_score`` (if given), or if the score has not improved over the last ``max_wait_count`` trials.
        
Verbosity levels:
    0: print nothing
    1: print only when next sequence is found
    2: print every step
'''
        if maximize:
            coef = 1
        else:
            coef = -1

        score = coef*self.scorer.score(seq)
    
        i = 0
        wait_count = 0
        while True:
            i += 1
            wait_count += 1 
            new_seq = self.change_random_codon(seq, start_window, end_window )
            new_score = coef*self.scorer.score(new_seq)
            if verbosity > 1:
                print(i, coef*score, coef*new_score)
            if new_score > score:
                score = new_score
                seq = new_seq
                wait_count = 0
                if verbosity > 0:
                    print(i, coef*score, seq)
                if target_score:
                    if abs(score-coef*target_score) <= tolerance:
                        break
            elif wait_count > max_wait_count:
                break
#        print(i, coef*score, seq)
        return (seq, score)
    
    
class TestScorer:
    def score(self, seq):
        i = 0
        for c in seq:
            if c == 'A':
                i += 1
        return i
    
def test():
    s = TestScorer()
    o = CodonOptimizer(s)

    seq = Seq('ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG', IUPAC.unambiguous_dna)

    o.hillclimb(seq, target_score=168, maximize=False)
        
# when run as its own script, 
if __name__ == "__main__":
    test()
