#!/usr/bin/env python3
import random, sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from codon_tools import *


class MultiScorer(CodonFreqScorer):
    def __init__(self, SS_coef, m2stop_coef, CpG_coef, UpA_coef):
        super().__init__()
        self.SS_coef = SS_coef
        self.m2stop_coef = m2stop_coef
        self.CpG_coef = CpG_coef
        self.UpA_coef = UpA_coef
        self.have_target_m2stop = False
        self.have_target_CpG = False
        self.have_target_UpA = False
        self.m2stop_count = 0
        self.CpG_count = 0
        self.UpA_count = 0
        
    def set_score_targets(self, m2stop_count = None, CpG_count = None, UpA_count = None):
        if m2stop_count is None:
            self.have_target_m2stop = False
        else:
            self.have_target_m2stop = True
            self.m2stop_count = m2stop_count
        
        if CpG_count is None:
            self.have_target_CpG = False
        else:
            self.have_target_CpG = True
            self.CpG_count = CpG_count
        
        if UpA_count is None:
            self.have_target_UpA = False
        else:
            self.have_target_UpA = True
            self.UpA_count = UpA_count
           
    def set_score_coefs(self, SS_coef = None, m2stop_coef = None,
                              CpG_coef = None, UpA_coef = None):
        if SS_coef is not None:
            self.SS_coef = SS_coef
        if m2stop_coef is not None:
            self.m2stop_coef = m2stop_coef
        if CpG_coef is not None:
            self.CpG_coef = CpG_coef
        if UpA_coef is not None:
            self.UpA_coef = UpA_coef


    def calc_score_components(self, seq):
        codon_freqs = self.sa.calc_syn_codon_freqs(seq)
        codon_SS = self.calc_SS_to_target(codon_freqs)
        m2stop_count = self.sa.count_muts_to_stop(seq)
        CpG_count, UpA_count = self.sa.count_CpG(seq)
        return codon_SS, m2stop_count, CpG_count, UpA_count

    def score(self, seq):
        codon_SS, m2stop_count, CpG_count, UpA_count = self.calc_score_components(seq)
        if self.have_target_m2stop:
            m2stop_score = (self.m2stop_count - m2stop_count)**2
        else:
            m2stop_score = m2stop_count
        if self.have_target_CpG:
            CpG_score = (self.CpG_count - CpG_count)**2
        else:
            CpG_score = CpG_count
        if self.have_target_UpA:
            UpA_score = (self.UpA_count - UpA_count)**2
        else:
            UpA_score = UpA_count
        return self.SS_coef*codon_SS + self.m2stop_coef*m2stop_score + self.CpG_coef*CpG_score + self.UpA_coef*UpA_score



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


def characterize_sequences(records):
    sa = SequenceAnalyzer()
    for record in records:
        print("Construct:", record.id)
        print("Description:", record.description)
        print("Length (# codons):", int(len(record.seq)/3))
        m2stop = sa.count_muts_to_stop(record.seq)
        (CpG, UpA) = sa.count_CpG(record.seq)
        print("Mutations to stop: %i, CpG count: %i, UpA count: %i" % (m2stop, CpG, UpA))
        print("Codon frequencies:")
        sa.calc_syn_codon_freqs(record.seq, verbosity = 1)
        print("-----------------------------------------------")


def generate_constructs(to_optimize, target):
    max_wait_count = 10000
    target_seq = target.seq
    seq = to_optimize.seq
    records = [to_optimize]
    
    records = records + \
                [match_codon_freqs(seq,
                  target_seq,
                  "C" + str(i),
                  "codon frequencies matched to Mahoney strain",
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
                  "codon frequencies matched to Mahoney strain and minimized mutations to stop, CpG, UpA",
                  max_wait_count) for i in [5, 6, 7]]

    records = records + \
                [match_all_freqs(seq,
                  target_seq,
                  "C" + str(i),
                  "codon frequencies and frequencies of mutations to stop, CpG, UpA all matched to Mahoney strain",
                  max_wait_count) for i in [8, 9, 10]]

    records = records + \
                [match_all_freqs_but_stop(seq,
                  target_seq,
                  "C" + str(i),
                  "codon frequencies and frequencies of CpG, UpA matched to Mahoney strain, minimized mutations to stop",
                  max_wait_count) for i in [11, 12, 13]]

    records = records + \
                [match_codon_freqs_minimize_CpG(seq,
                  target_seq,
                  "C" + str(i),
                  "codon frequencies matched to Mahoney strain, minimized CpG only",
                  max_wait_count) for i in [14, 15, 16]]


    characterize_sequences([target] + records)
    with open("out.fasta", 'w') as outfile:
        SeqIO.write(records, outfile, "fasta")
   
   
   
# when run as its own script
if __name__ == "__main__":
    random.seed(123)
    # read in the sequence defining target codon frequencies    
    target = SeqIO.parse(open("fasta/Mahoney_ORF.fasta", "rU"), "fasta").__next__()

    # read in the sequence to optimize
    to_optimize = SeqIO.parse(open("fasta/GFP_wt.fasta", "rU"), "fasta").__next__()

    # do the job  
    generate_constructs(to_optimize, target)
