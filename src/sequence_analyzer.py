from codon_analyzer import CodonAnalyzer

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
        

def test():
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC

    sa = SequenceAnalyzer()
    seq = Seq('ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG', IUPAC.unambiguous_dna)

    m2stop_count = sa.count_muts_to_stop(seq)
    print(m2stop_count)
    CpG_count, UpA_count = sa.count_CpG(seq)
    print(CpG_count, UpA_count)

    
# when run as its own script, 
if __name__ == "__main__":
    test()
