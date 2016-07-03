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


    def print_codon_freqs(self, codon_freqs, ndigits = 3):
        for aa in reverse_genetic_code:
            if aa == '*': # we skip over stop codons
                continue
            print(aa + ": ", end='')
            codons = reverse_genetic_code[aa]
            fam_freqs = {c:codon_freqs[c] for c in codons}
            for codon in fam_freqs:
                print(codon + ":", round(fam_freqs[codon], 3), end='; ')
            print()


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

    def calc_syn_codon_freqs(self, seq, verbosity = 0):
        """Calculates the relative frequencies of different codons in the same codon family. Stop codons are ignored. Returns a dictionary with relative frequencies, normalized such that they sum to one within each amino-acid family (i.e., the total sum over all families is one).
"""
        counts = self.count_codons(seq) # get absolute codon counts
        if verbosity > 1:
            print('--Absolute counts--')
            self.print_codon_freqs(counts)

        
        # normalize by amino-acid family
        for aa in reverse_genetic_code:
            if aa == '*': # we skip over stop codons
                continue
            codons = reverse_genetic_code[aa]
            fam_counts = {c:counts[c] for c in codons}
            total = sum(fam_counts.values())
            for c in fam_counts:
                counts[c] = fam_counts[c]/total
        
        if verbosity > 0:
            print('--Relative frequencies--')
            self.print_codon_freqs(counts)


def test(seq):
    sa = SA()
    sa.calc_syn_codon_freqs(seq, verbosity = 2)
    
    
        
# when run as its own script, 
if __name__ == "__main__":
    # read in the sequence    
    seq = SeqIO.parse(open("fasta/Mahoney_ORF.fasta", "rU"), "fasta").__next__().seq
    # run analysis    
    test(seq)
    print()
    
    # read in the sequence    
    seq = SeqIO.parse(open("fasta/GFP_wt.fasta", "rU"), "fasta").__next__().seq
    # run analysis    
    test(seq)
