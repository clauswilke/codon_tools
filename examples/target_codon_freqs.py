#!/usr/bin/env python3
from codon_tools import SequenceAnalyzer, CodonOptimizer, StopAndCpGScorer
from codon_tools import reverse_genetic_code

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


class SA(SequenceAnalyzer):
    def __init__(self):
        self.zero_codon_counts = {}
        for c1 in ['A', 'G', 'C', 'T']:
            for c2 in ['A', 'G', 'C', 'T']:
                for c3 in ['A', 'G', 'C', 'T']:
                    codon = c1 + c2 + c3
                    # no stop codons
                    if codon in ['TAA', 'TAG', 'TGA']:
                        continue
                    self.zero_codon_counts[codon] = 0



    def count_codons(self, seq):
        """Counts how many times each codon appears in a sequence, and returns a dictionary with codon counts.
"""
        assert len(seq) % 3 == 0
        
        counts = self.zero_codon_counts.copy()
        for i in range(int(len(seq)/3)):
            codon = str(seq[3*i:3*i+3])
            if codon in counts:
                counts[codon] += 1
            # we ignore any codons that we cannot uniquely identify
        return counts

    def calc_syn_codon_freqs(self, seq):
        """Calculates the relative frequencies of different codons in the same codon family. Stop codons are ignored. Returns a dictionary with relative frequencies, normalized such that they sum to one within each amino-acid family (i.e., the total sum over all families is one).
"""
        counts = self.count_codons(seq) # get absolute codon counts
        print(counts)
        # normalize by amino-acid family
        for aa in reverse_genetic_code:
            if aa == '*': # we skip over stop codons
                continue
            codons = reverse_genetic_code[aa]
            fam_counts = {c:counts[c] for c in codons}
            total = sum(fam_counts.values())
            for c in fam_counts:
                counts[c] = fam_counts[c]/total
        print(counts)
        print(sum(counts.values()))
        for aa in reverse_genetic_code:
            if aa == '*': # we skip over stop codons
                continue
            print(aa)
            codons = reverse_genetic_code[aa]
            fam_counts = {c:counts[c] for c in codons}
            print(fam_counts)
            print()


def test(seq):
    sa = SA()
    sa.calc_syn_codon_freqs(seq)
    

def de_CpG(seq):
    """Remove as many CpGs as possible, as well as potential single-point mutations to stop and UpAs.
"""
    scorer = StopAndCpGScorer(1, 1, 1)
    o = CodonOptimizer(scorer)

    seq_orig = seq
    seq, score = o.hillclimb(seq, maximize=False, max_wait_count = 5000, verbosity = 0)
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
    seq = SeqIO.parse(open("fasta/GFP_wt.fasta", "rU"), "fasta").__next__().seq
    # run analysis    
    test(seq)
