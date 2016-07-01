#!/usr/bin/env python3
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord    
    
from codon_tools import CodonOptimizer, FopScorer, opt_codons_E_coli


def deoptimize(seq, gene_description, start_window = 14, end_window = 14, opt_codons = opt_codons_E_coli):
    """Parameters:
* ``seq``: coding sequence to de-optimize
* ``start_window``: number of codons to omit at beginning of sequence
* ``end_window``: number of codons to omit at end of sequence
* ``opt_codons``: set of optimal codons

"""
    assert len(seq) % 3 == 0
    
    scorer = FopScorer()
    o = CodonOptimizer(scorer)

    seq_orig = seq
    tolerance = 0.01 # how closely do we want to match the final number?
    for i in range(50, 5, -10):
        seq, score = o.hillclimb(seq, start_window, end_window,
                        i/100-tolerance, tolerance, maximize=False,
                        max_wait_count = 5000, verbosity = 0)
        assert seq_orig.translate() == seq.translate()

        description = gene_description + " -- deoptimized to Fop = %f (keeping first %i and last %i codons unchanged)" % (score, start_window, end_window) 
        seq_record = SeqRecord( seq, id='', description=description )
        print(seq_record.format("fasta"))

    
        
# when run as its own script, 
if __name__ == "__main__":
    # set up command line arguments
    parser = argparse.ArgumentParser(description='Codon-deoptimize gene sequence.')
    parser.add_argument('filename', metavar='filename',
                        help='fasta-formatted file holding the sequence')
    args = parser.parse_args()

    # read in the sequence     
    record = SeqIO.parse(open(args.filename, "rU"), "fasta").__next__()
    # run analysis    
    deoptimize(record.seq, record.description, 14, 14)
