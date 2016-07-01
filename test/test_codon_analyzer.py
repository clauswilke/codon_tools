#!/usr/bin/env python3
import time       

from codon_tools import CodonAnalyzer

def test_stop_count():
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
