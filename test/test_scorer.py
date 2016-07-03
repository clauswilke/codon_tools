def test(seq, target_seq):
    scorer = CodonFreqScorer()
    scorer.set_codon_freqs_from_seq(target_seq)
    print(scorer.score(seq))
