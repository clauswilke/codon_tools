import random

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna

class CodonShuffler:
    """
    Class to shuffle codons within an ORF, without changing the overall
    codon frequency.
    """
    def __init__(self):
        pass

    def make_lookup_table(self, seq):
        """
        Make a reverse codon lookup table that lists codons used in seq for
        each amino acid. The same codon may be listed more than once for each
        amino acid.
        """
        codons = {}
        seq = Seq(seq, generic_dna)

        for i in range(int(len(seq)/3)):
            codon = seq[3*i:3*i+3]
            amino_acid = str(Seq.translate(codon))
            if amino_acid in codons:
                codons[amino_acid].append(str(codon))
            else:
                codons[amino_acid] = [str(codon)]

        return codons

    def shuffle_codons(self, seq):
        """
        Shuffle codons in a sequence without changing codon frequencies.
        """
        # Construct dictionary of all codons present in the sequence
        codons = self.make_lookup_table(seq)
        # Shuffle codons
        for amino_acid in codons:
            random.shuffle(codons[amino_acid])
        # Translate original sequence
        seq_aa = Seq(seq, generic_dna).translate()
        shuffled_seq = ""
        # Reconstruct sequence with codons shuffled
        for amino_acid in str(seq_aa):
            shuffled_seq += codons[amino_acid].pop()

        return shuffled_seq
