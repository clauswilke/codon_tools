# Optimal codons for E. coli
# T. Zhou, M. Weems, C. O. Wilke (2009). Translationally optimal codons associate with structurally sensitive sites in proteins. Mol. Biol. Evol. 26:1571–1580. 
opt_codons_E_coli = { 'A':['GCT'], 'R':['CGT', 'CGC'], 'N':['AAC'], 'D':['GAC'], 'C':['TGC'], 'Q':['CAG'], 'E':['GAA'], 'G':['GGT','GGC'], 'H':['CAC'], 'I':['ATC'], 'L':['CTG'], 'F':['TTC'], 'P':['CCG'], 'S':['TCT','TCC'], 'T':['ACT','ACC'], 'Y':['TAC'], 'V':['GTT','GTA'] }


# Reverse genetic code
# For each amino acid, provides a list of all the codons that translate into
# that amino acid
reverse_genetic_code = {  'A':['GCA', 'GCC', 'GCG', 'GCT'], 'R':['AGA', 'AGG', 'CGA', 'CGT', 'CGC', 'CGG'], 'N':['AAC', 'AAT'], 'D':['GAC', 'GAT'], 'C':['TGC', 'TGT' ], 'Q':['CAA', 'CAG'], 'E':['GAA', 'GAG'], 'G':['GGA', 'GGC', 'GGG', 'GGT'], 'H':['CAC', 'CAT'], 'I':['ATA', 'ATC', 'ATT'], 'L':['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'], 'F':['TTT', 'TTC'], 'P':['CCA', 'CCC', 'CCG', 'CCT'], 'S':['AGC', 'AGT', 'TCA', 'TCC','TCG','TCT'], 'T':[ 'ACA', 'ACT', 'ACC', 'ACG' ], 'Y':['TAC', 'TAT'], 'V':['GTA', 'GTC', 'GTG', 'GTT'], 'W':['TGG'], 'M':['ATG'], 'K':['AAA', 'AAG'], '*':['TAA', 'TAG', 'TGA']}

