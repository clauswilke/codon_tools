#!/Users/wilke/anaconda/bin/python3
from sequence_analyzer import SequenceAnalyzer
from codon_analyzer import CodonAnalyzer
from codon_optimizer import CodonOptimizer, reverse_genetic_code

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

GFP_wt_seq = Seq('ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG', IUPAC.unambiguous_dna)


class StopAndCpGScorer:
    def __init__(self, m2stop_coef, CpG_coef, UpA_coef):
        self.sa = SequenceAnalyzer()
        self.m2stop_coef = m2stop_coef
        self.CpG_coef = CpG_coef
        self.UpA_coef = UpA_coef

    def calc_score_components(self, seq):
        m2stop_count = self.sa.count_muts_to_stop(seq)
        CpG_count, UpA_count = self.sa.count_CpG(seq)
        return m2stop_count, CpG_count, UpA_count

    def score(self, seq):
        m2stop_count, CpG_count, UpA_count = self.calc_score_components(seq)
        return self.m2stop_coef*m2stop_count + self.CpG_coef*CpG_count + self.UpA_coef*UpA_count
    
    
def analyze_m2stop():
    ca = CodonAnalyzer()
    
    seq = GFP_wt_seq
    
    print("""Each line contains one codon of wt GFP, the corresponding amino acid,
and the number of single-point mutations that lead to a stop codon.
For codons for which the count exceeds 0, all possible alternative codons
and their counts are shown as well.""")
    for i in range(int(len(seq)/3)):
        codon = seq[3*i:3*i+3]
        count = ca.count_stop_muts(codon)
        aa = str(codon.translate())
        print(codon, aa, count)
        if count > 0:
            all_codons = reverse_genetic_code[aa]
            print("alternative codons for " + aa + ":")
            for c in all_codons:
                print(" ", c, ca.count_stop_muts(c))
        
def main():
    scorer = StopAndCpGScorer(1, 1, 1)
    o = CodonOptimizer(scorer)

    seq = GFP_wt_seq
    seq_orig = seq
    seq, score = o.hillclimb(seq, maximize=False, max_wait_count = 5000, verbosity=0)
    assert seq_orig.translate() == seq.translate()
    
    m2stop_count, CpG_count, UpA_count = scorer.calc_score_components(seq_orig)
    print("Original sequence")
    print("=================")
    print("Mutations to stop: %i, CpG count: %i, UpA count: %i" % (m2stop_count, CpG_count, UpA_count))

    m2stop_count, CpG_count, UpA_count = scorer.calc_score_components(seq)    
    print("\nOptimized sequence:")
    print("=================")
    print("Mutations to stop: %i, CpG count: %i, UpA count: %i\n" % (m2stop_count, CpG_count, UpA_count))
    print(seq)
    
        
# when run as its own script, 
if __name__ == "__main__":
    print(len(GFP_wt_seq))
    #main()
    #analyze_m2stop()
  
# Mahoney CpGs: 174/6627, 2.6%
#         UpAs: 363/6627, 5.5%
  
# GFP CpGs: 60/717, 8.4%
#     UpAs: 13/717, 1.8%
