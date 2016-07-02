#!/usr/bin/env python3
import argparse, sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord    
    
from codon_tools import CodonOptimizer, FopScorer, opt_codons_E_coli


def deoptimize(seq, gene_description, Fop_start, Fop_step, Fop_stop, start_window = 14, end_window = 14, max_wait_count = 5000, opt_codons = opt_codons_E_coli):
    """Parameters:
* ``seq``: coding sequence to de-optimize
* ``start_window``: number of codons to omit at beginning of sequence
* ``end_window``: number of codons to omit at end of sequence
* ``opt_codons``: set of optimal codons

"""
    assert len(seq) % 3 == 0
    
    scorer = FopScorer()
    o = CodonOptimizer(scorer)

    score = scorer.score(seq)
    seq_orig = seq
    Fop_orig = score
    tolerance = 0.01 # how closely do we want to match the final number?
    records = [SeqRecord(seq, id='',
                 description = gene_description + " -- Fop = %f" % Fop_orig)]
    while True:
        if Fop_step < 0:
            Fop_target = Fop_start - tolerance
        else:
            Fop_target = Fop_start + tolerance
        if Fop_target < score:
            maximize = False
        else:
            maximize = True
        seq, score = o.hillclimb(seq, start_window, end_window,
                        Fop_target, tolerance, max_wait_count,
                        maximize, verbosity = 0)
        assert seq_orig.translate() == seq.translate()
        
        description = gene_description + " -- recoded to Fop = %f (keeping first %i and last %i codons unchanged)" % (score, start_window, end_window) 
        seq_record = SeqRecord(seq, id = '', description = description)
        records.append(seq_record)
        
        Fop_start += Fop_step
        if Fop_step*(Fop_stop-Fop_start) < 0:
            break
    
    return records
    
        
# when run as its own script, 
if __name__ == "__main__":
    # set up command line arguments
    parser = argparse.ArgumentParser(description='Codon-deoptimize gene sequence.')
    parser.add_argument('filename', metavar='filename',
                        help='fasta-formatted file holding the sequence')
    parser.add_argument('-x', '--exclude-front', default=14, type=int,
                        metavar='n',
                        help='number of codons to exclude from the front of the sequence, default = 14')
    parser.add_argument('-X', '--exclude-back', default=14, type=int,
                        metavar='n',
                        help='number of codons to exclude from the back of the sequence, default = 14')
    parser.add_argument('-a', '--Fop-start', default=0.5, type=float,
                        metavar='x',
                        help='starting value for Fop, default = 0.5')
    parser.add_argument('-d', '--Fop-step', default=-0.1, type=float,
                        metavar='x',
                        help='step size from one Fop value to the next, default = -0.1')
    parser.add_argument('-e', '--Fop-stop', default=0.1, type=float,
                        metavar='x',
                        help='end value for Fop, default = 0.1')

    parser.add_argument('-M', '--max-wait', default=5000, type=int,
                        metavar='n',
                        help='maximum number of attempted codon substitutions before we give up the optimization procedure, default = 5000')


    args = parser.parse_args()

    # read in the sequence     
    record = SeqIO.parse(open(args.filename, "rU"), "fasta").__next__()
    # run analysis    
    records = deoptimize(record.seq, record.description,
                         args.Fop_start, args.Fop_step, args.Fop_stop,
                         args.exclude_front, args.exclude_back,
                         args.max_wait)
    SeqIO.write(records, sys.stdout, "fasta")
