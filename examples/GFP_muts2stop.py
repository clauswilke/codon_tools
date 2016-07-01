#!/usr/bin/env python3
from codon_tools import CodonAnalyzer, reverse_genetic_code

from Bio import SeqIO
from Bio.Seq import Seq
    
def analyze_m2stop(seq):
    print("""Each line contains one codon of wt GFP, the corresponding amino acid,
and the number of single-point mutations that lead to a stop codon.
For codons for which the count exceeds 0, all possible alternative codons
and their counts are shown as well.""")

    ca = CodonAnalyzer()

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
            
        
# when run as its own script, 
if __name__ == "__main__":
    # read in GFP wt sequence    
    seq = SeqIO.parse(open("GFP_wt.fasta", "rU"), "fasta").__next__().seq
    # run analysis
    analyze_m2stop(seq)
