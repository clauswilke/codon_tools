from .sequence_analyzer import SequenceAnalyzer
from .lookup_tables import opt_codons_E_coli

class StopAndCpGScorer:
    """Scores sequences based on the number of possible mutations to stop codons, number of CpGs, and number of UpAs. In the constructor, it is possible to set the relative weight of these three counts, including setting some of them to zero.
"""
    def __init__(self, m2stop_coef, CpG_coef, UpA_coef):
        self.sa = SequenceAnalyzer()
        self.m2stop_coef = m2stop_coef
        self.CpG_coef = CpG_coef
        self.UpA_coef = UpA_coef

    def calc_score_components(self, seq):
        m2stop_count = self.sa.count_muts_to_stop(seq)
        CpG_count, UpA_count = self.sa.count_CpG(seq)
        return m2stop_count, CpG_count, UpA_count

    def score(self, seq):
        m2stop_count, CpG_count, UpA_count = self.calc_score_components(seq)
        return self.m2stop_coef*m2stop_count + self.CpG_coef*CpG_count + self.UpA_coef*UpA_count


class FopScorer:
    """Scores sequences based on the fraction of optimal codons. By default, uses the E. coli optimal codons from Zhou et al. 2009. A different set of optimal codons can be specified in the constructor.
"""
    def __init__(self, opt_codons = opt_codons_E_coli):
        self.sa = SequenceAnalyzer()
        self.opt_codons = opt_codons

    def score(self, seq):
        opt_count, total_count, Fop = self.sa.calc_Fop(seq, self.opt_codons)
        return Fop


class CodonFreqScorer:
    """Scores sequences based on how close codon frequencies are to target frequencies.
"""
    def __init__(self):
        self.have_codon_freqs = False
        self.sa = SequenceAnalyzer()

    def set_codon_freqs_from_seq(self, seq):
        """Sets target codon frequencies calculated from the given sequence.
"""
        self.codon_freqs = self.sa.calc_syn_codon_freqs(seq)
        self.have_codon_freqs = True 
        
    def calc_SS_to_target(self, codon_freqs):
        """Calculates the sum of the square deviations between the provided codon frequencies and the target frequencies. Returns 0 if target frequencies have not been set.
"""
        if not self.have_codon_freqs:
            return 0
        
        sum_square = 0
        for codon in codon_freqs:
            sum_square += (self.codon_freqs[codon] - codon_freqs[codon])**2
        return sum_square

    def score(self, seq):
        codon_freqs = self.sa.calc_syn_codon_freqs(seq)
        return(self.calc_SS_to_target(codon_freqs))



class MultiScorer(CodonFreqScorer):
    """Builds on the ``CodonFreqScorer`` but can also score the number of point mutations to a stop codon, CpGs, and UpAs. For the latter three quantities, can score towards a target value or towards minimization/maximization.
"""
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
