#! /usr/bin/env python3
import sys

from .codon_analyzer import *
from .lookup_tables import opt_codons_E_coli, reverse_genetic_code

class SequenceAnalyzer:
    def __init__(self):
        self.codon_analyzer = CodonAnalyzer()

        self.zero_codon_counts = {}
        for c1 in ['A', 'G', 'C', 'T']:
            for c2 in ['A', 'G', 'C', 'T']:
                for c3 in ['A', 'G', 'C', 'T']:
                    codon = c1 + c2 + c3
                    # no stop codons
                    if codon in ['TAA', 'TAG', 'TGA']:
                        continue
                    self.zero_codon_counts[codon] = 0



    def count_CpG(self, seq):
        """Counts the number of CpGs and UpAs in a sequence.
Note: The code expects Ts not Us in the sequence.
"""
        CpG_count = 0
        UpA_count = 0
        for i in range(len(seq)-1):
            if seq[i] == 'C':
                if seq[i+1] == 'G':
                    CpG_count += 1
            if seq[i] == 'T':
                if seq[i+1] == 'A':
                    UpA_count += 1
        return (CpG_count, UpA_count)

    def count_muts_to_stop(self, seq):
        """Counts the number of possible single-point mutations to a stop codon.
"""
        assert len(seq) % 3 == 0

        count = 0
        for i in range(int(len(seq)/3)):
            codon = seq[3*i:3*i+3]
            count += self.codon_analyzer.count_stop_muts(codon)
        return count

    def calc_Fop(self, seq, opt_codons = opt_codons_E_coli, verbosity = 0):
        """Calculates the fraction of optimal codons in a sequence. The fraction is calculated relative to the number of sites where an optimal codon is possible (i.e., for example, methionines are excluded since they are coded by only one amino acid).
"""
        assert len(seq) % 3 == 0

        codon_list = [j for i in opt_codons.values() for j in i]
        aa_list = opt_codons.keys()
        total_count = 0
        opt_count = 0

        for i in range(int(len(seq)/3)):
            codon = seq[3*i:3*i+3]
            aa = str( codon.translate() )
            ignore = True
            opt = False
            if aa in aa_list:
                ignore = False
                c = str( codon )
                opt = (c in codon_list)
            total_count += not ignore
            opt_count += opt

            if verbosity > 1:
                print(3*i, aa, str( codon ), end=' ')
                if ignore:
                    print("excluded")
                elif opt:
                    print("optimal")
                else:
                    print()
        if verbosity > 1:
            print("\nSummary:")
        if verbosity > 0:
            print("Count of optimal codons =", opt_count)
            print("Count of non-excluded sites =", total_count)
            print("Fraction of optimal codons =", float(opt_count)/total_count)
        return ( opt_count, total_count, float(opt_count)/total_count )


    def print_codon_freqs(self, codon_freqs, outfile = sys.stdout, ndigits = 3):
        for aa in reverse_genetic_code:
            if aa == '*': # we skip over stop codons
                continue
            outfile.write(aa + ": ")
            codons = reverse_genetic_code[aa]
            fam_freqs = {c:codon_freqs[c] for c in codons}
            for codon in fam_freqs:
                outfile.write("%s: %f" % (codon, round(fam_freqs[codon], 3)))
                outfile.write(";")
            outfile.write("\n")


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

    def calc_syn_codon_freqs(self, seq, outfile = sys.stdout, verbosity = 0):
        """Calculates the relative frequencies of different codons in the same codon family. Stop codons are ignored. Returns a dictionary with relative frequencies, normalized such that they sum to one within each amino-acid family (i.e., the total sum over all families is one).
"""
        counts = self.count_codons(seq) # get absolute codon counts
        if verbosity > 1:
            outfile.write('--Absolute counts--\n')
            self.print_codon_freqs(counts, outfile)


        # normalize by amino-acid family
        for aa in reverse_genetic_code:
            if aa == '*': # we skip over stop codons
                continue
            codons = reverse_genetic_code[aa]
            fam_counts = {c:counts[c] for c in codons}
            total = sum(fam_counts.values())
            if total > 0:
                for c in fam_counts:
                    counts[c] = fam_counts[c]/total

        if verbosity > 0:
            if verbosity > 1:
                outfile.write('--Relative frequencies--\n')
            self.print_codon_freqs(counts, outfile)
        return counts

    def count_diffs(self, seq1, seq2):
        """Count the number of differences between two sequences of the same
        length."""
        assert len(seq1) == len(seq2)
        diffs = sum(1 for a, b in zip(seq1, seq2) if a != b)

        return diffs
