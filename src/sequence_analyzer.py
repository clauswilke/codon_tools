#! /usr/bin/env python3
from .codon_analyzer import *
from .lookup_tables import opt_codons_E_coli

class SequenceAnalyzer:
    def __init__(self):
        self.codon_analyzer = CodonAnalyzer()

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
