#!/usr/bin/env python3
from codon_tools import SequenceAnalyzer, CodonOptimizer

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

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
    
        
def de_CpG(seq):
    """Remove as many CpGs as possible, as well as potential single-point mutations to stop and UpAs.
"""
    scorer = StopAndCpGScorer(1, 1, 1)
    o = CodonOptimizer(scorer)

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
    # read in GFP wt sequence    
    seq = SeqIO.parse(open("GFP_wt.fasta", "rU"), "fasta").__next__().seq
    # run analysis    
    de_CpG(seq)
