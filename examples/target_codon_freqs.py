#!/usr/bin/env python3
from codon_tools import SequenceAnalyzer, CodonOptimizer, CodonFreqScorer, StopAndCpGScorer

from Bio import SeqIO

def optimize_codon_freqs(seq, target_seq):
    scorer = CodonFreqScorer()
    scorer.set_codon_freqs_from_seq(target_seq)
    o = CodonOptimizer(scorer)
    sa = SequenceAnalyzer()

    seq_orig = seq
    seq, score = o.hillclimb(seq, maximize=False, max_wait_count = 5000, verbosity = 0)
    assert seq_orig.translate() == seq.translate()
    
    print("Original sequence")
    print("=================")
    sa.calc_syn_codon_freqs(seq_orig, verbosity = 1)

    print("\nTarget sequence")
    print("=================")
    sa.calc_syn_codon_freqs(target_seq, verbosity = 1)

    print("\nOptimized sequence:")
    print("=================")
    sa.calc_syn_codon_freqs(seq, verbosity = 1)
    
    print()
    print(seq)
    
    scorer2 = StopAndCpGScorer(1, 1, 1)
    m2stop_count, CpG_count, UpA_count = scorer2.calc_score_components(seq)
    print("Mutations to stop: %i, CpG count: %i, UpA count: %i" % (m2stop_count, CpG_count, UpA_count))

 
        
# when run as its own script, 
if __name__ == "__main__":
    # read in the sequence defining target codon frequencies    
    target_seq = SeqIO.parse(open("fasta/Mahoney_ORF.fasta", "rU"), "fasta").__next__().seq

    # read in the sequence to optimize
    seq = SeqIO.parse(open("fasta/GFP_wt.fasta", "rU"), "fasta").__next__().seq

    # do the job  
    optimize_codon_freqs(seq, target_seq)
