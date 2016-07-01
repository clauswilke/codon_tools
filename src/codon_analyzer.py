#! /usr/bin/env python3
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

class CodonAnalyzer:
    def __init__(self):
        self.stop_counts = {} # initialize stop_counts dict as empty, it will be filled as needed


    def count_stop_muts(self, codon):
        '''Function that takes a single codon as argument (string or Bio Seq object) and returns the number of single-point mutations that lead to stop codons.'''
        codon = str(codon)  # make sure we work with strings, not biopython sequences
        assert len(codon) == 3
        
        if codon in self.stop_counts:
            return self.stop_counts[codon]
        
        stop_count = 0
        for c in ['A', 'G', 'C', 'T']:
            # position 1
            if codon[0] != c:
                codon_mut = c + codon[1:3]
#                print(codon, "->", codon_mut)
                if Seq(codon_mut, generic_dna).translate() == '*':
                    stop_count += 1
#                    print("  mutation to stop")

            if codon[1] != c:
                codon_mut = codon[0] + c + codon[2]
#                print(codon, "->", codon_mut)
                if Seq(codon_mut, generic_dna).translate() == '*':
                    stop_count += 1
#                    print("  mutation to stop")

            if codon[2] != c:
                codon_mut = codon[0:2] + c
#                print(codon, "->", codon_mut)
                if Seq(codon_mut, generic_dna).translate() == '*':
                    stop_count += 1
#                    print("  mutation to stop")

        self.stop_counts[codon] = stop_count # save for future use
        return stop_count       
        

def test_stop_count():
    import time       
    a = CodonAnalyzer()

    t0 = time.time()
    for c1 in ['A', 'G', 'C', 'T']:
        for c2 in ['A', 'G', 'C', 'T']:
            for c3 in ['A', 'G', 'C', 'T']:
                codon = c1+c2+c3
                count = a.count_stop_muts(codon)
                print(codon + ": " + str(count))
            
    t1 = time.time()
    for c1 in ['A', 'G', 'C', 'T']:
        for c2 in ['A', 'G', 'C', 'T']:
            for c3 in ['A', 'G', 'C', 'T']:
                codon = c1+c2+c3
                count = a.count_stop_muts(codon)
                print(codon + ": " + str(count))

    t2 = time.time()
    print(t1-t0, t2-t1)
    
# when run as its own script, 
if __name__ == "__main__":
    test_stop_count()
