#!/usr/bin/env python3
from codon_tools import SequenceAnalyzer

from Bio import SeqIO


def test(seq):
    sa = SequenceAnalyzer()
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
