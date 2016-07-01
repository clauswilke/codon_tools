#! /usr/bin/env python3
from .codon_analyzer import *

class SequenceAnalyzer:
    def __init__(self):
        self.codon_analyzer = CodonAnalyzer()

    def count_CpG(self, seq):
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
        assert len(seq) % 3 == 0
        
        count = 0
        for i in range(int(len(seq)/3)):
            codon = seq[3*i:3*i+3]
            count += self.codon_analyzer.count_stop_muts(codon)
        return count
        
