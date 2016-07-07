#!/usr/bin/env python3
import random, sys, argparse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from codon_tools import *



def match_codon_freqs(seq, target_seq, id, description, max_wait_count):
    scorer = CodonFreqScorer()
    scorer.set_codon_freqs_from_seq(target_seq)
    o = CodonOptimizer(scorer)
    
    seq_orig = seq
    seq, score = o.hillclimb(seq, maximize=False, max_wait_count = max_wait_count, verbosity = 0)
    assert seq_orig.translate() == seq.translate()
    return SeqRecord(seq, id = id, description = description)    


def minimize_stop_CpG(seq, id, description, max_wait_count):
    scorer = StopAndCpGScorer(1, 1, 1)
    o = CodonOptimizer(scorer)
    
    seq_orig = seq
    seq, score = o.hillclimb(seq, maximize=False, max_wait_count = max_wait_count, verbosity = 0)
    assert seq_orig.translate() == seq.translate()
    return SeqRecord(seq, id = id, description = description)    


def match_codon_freqs_minimize_stop_CpG(seq, target_seq, id, description, max_wait_count):
    scorer = MultiScorer(1, .1, .1, .1)
    scorer.set_codon_freqs_from_seq(target_seq)
    o = CodonOptimizer(scorer)
    
    seq_orig = seq
    seq, score = o.hillclimb(seq, maximize=False, max_wait_count = max_wait_count, verbosity = 0)
    assert seq_orig.translate() == seq.translate()
    return SeqRecord(seq, id = id, description = description)    


def match_all_freqs(seq, target_seq, id, description, max_wait_count):
    scorer = MultiScorer(1, 1, 1, 1)
    scorer.set_codon_freqs_from_seq(target_seq)
    codon_SS, m2stop_count, CpG_count, UpA_count = scorer.calc_score_components(target_seq)
    scorer.set_score_targets(m2stop_count = m2stop_count*len(seq)/len(target_seq),
                             CpG_count = CpG_count*len(seq)/len(target_seq),
                             UpA_count = UpA_count*len(seq)/len(target_seq))    
    
    o = CodonOptimizer(scorer)
    
    seq_orig = seq
    seq, score = o.hillclimb(seq, maximize=False, max_wait_count = max_wait_count, verbosity = 0)
    assert seq_orig.translate() == seq.translate()
    return SeqRecord(seq, id = id, description = description)    


def match_codon_freqs_minimize_CpG(seq, target_seq, id, description, max_wait_count):
    scorer = MultiScorer(1, 0, 1, 0)
    scorer.set_codon_freqs_from_seq(target_seq)
    o = CodonOptimizer(scorer)
    
    seq_orig = seq
    seq, score = o.hillclimb(seq, maximize=False, max_wait_count = max_wait_count, verbosity = 0)
    assert seq_orig.translate() == seq.translate()
    return SeqRecord(seq, id = id, description = description)    



def match_all_freqs_but_stop(seq, target_seq, id, description, max_wait_count):
    scorer = MultiScorer(1, 1, 1, 1)
    scorer.set_codon_freqs_from_seq(target_seq)
    codon_SS, m2stop_count, CpG_count, UpA_count = scorer.calc_score_components(target_seq)
    scorer.set_score_targets(CpG_count = CpG_count*len(seq)/len(target_seq),
                             UpA_count = UpA_count*len(seq)/len(target_seq))    
    
    o = CodonOptimizer(scorer)
    
    seq_orig = seq
    seq, score = o.hillclimb(seq, maximize=False, max_wait_count = max_wait_count, verbosity = 0)
    assert seq_orig.translate() == seq.translate()
    return SeqRecord(seq, id = id, description = description)    


def characterize_sequences(records, outfile):
    sa = SequenceAnalyzer()
    for record in records:
        outfile.write("Construct: %s\n" % record.id)
        outfile.write("Description: %s\n" % record.description)
        outfile.write("Length (# codons): %i\n" % int(len(record.seq)/3))
        m2stop = sa.count_muts_to_stop(record.seq)
        (CpG, UpA) = sa.count_CpG(record.seq)
        outfile.write("Mutations to stop: %i, CpG count: %i, UpA count: %i\n" % (m2stop, CpG, UpA))
        outfile.write("Codon frequencies:\n")
        sa.calc_syn_codon_freqs(record.seq, outfile, verbosity = 2)
        outfile.write("-----------------------------------------------\n")


def generate_constructs(to_optimize, target, max_wait_count = 10000):
    target_seq = target.seq
    seq = to_optimize.seq
    records = [to_optimize]
    
    records = records + \
                [match_codon_freqs(seq,
                  target_seq,
                  "C" + str(i),
                  "codon frequencies matched to " + target.id + " strain",
                  max_wait_count) for i in [1, 2, 3]]

    records = records + \
                [minimize_stop_CpG(seq,
                  "C" + str(i),
                  "minimized mutations to stop, CpG, UpA",
                  max_wait_count) for i in [4]]

    records = records + \
                [match_codon_freqs_minimize_stop_CpG(seq,
                  target_seq,
                  "C" + str(i),
                  "codon frequencies matched to " + target.id + " strain and minimized mutations to stop, CpG, UpA",
                  max_wait_count) for i in [5, 6, 7]]

    records = records + \
                [match_all_freqs(seq,
                  target_seq,
                  "C" + str(i),
                  "codon frequencies and frequencies of mutations to stop, CpG, UpA all matched to " + target.id,
                  max_wait_count) for i in [8, 9, 10]]

    records = records + \
                [match_all_freqs_but_stop(seq,
                  target_seq,
                  "C" + str(i),
                  "codon frequencies and frequencies of CpG, UpA matched to " + target.id + ", minimized mutations to stop",
                  max_wait_count) for i in [11, 12, 13]]

    records = records + \
                [match_codon_freqs_minimize_CpG(seq,
                  target_seq,
                  "C" + str(i),
                  "codon frequencies matched to " + target.id + ", minimized CpG only",
                  max_wait_count) for i in [14, 15, 16]]

    fasta_fname = "recoded=%s_ref=%s.fasta" % (to_optimize.id, target.id)
    descr_fname = "recoded=%s_ref=%s_descr.txt" % (to_optimize.id, target.id)

    with open(descr_fname, 'w') as outfile:
        characterize_sequences([target] + records, outfile)
    
    with open(fasta_fname, 'w') as outfile:
        SeqIO.write(records, outfile, "fasta")
   
   
    
    
if __name__ == "__main__":
    # set up command line arguments
    parser = argparse.ArgumentParser(description='Recode gene sequence to reference')
    parser.add_argument('to_recode', metavar='to-recode',
                        help='fasta-formatted file holding the sequence to recode')
                        
    parser.add_argument('reference', metavar='reference',
                        help='fasta-formatted file holding the reference sequence')

    parser.add_argument('-M', '--max-wait', default=10000, type=int,
                        metavar='n',
                        help='maximum number of attempted codon substitutions before we give up the optimization procedure, default = 10000')

    parser.add_argument('-s', '--seed', default=123, type=int,
                        metavar='n',
                        help='random seed, default = 123')


    args = parser.parse_args()

    random.seed(args.seed)
    # read in the sequence defining target codon frequencies    
    target = SeqIO.parse(open(args.reference, "rU"), "fasta").__next__()

    # read in the sequence to optimize
    to_optimize = SeqIO.parse(open(args.to_recode, "rU"), "fasta").__next__()

    # do the job  
    generate_constructs(to_optimize, target, args.max_wait)
