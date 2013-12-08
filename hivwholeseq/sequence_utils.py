# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/12/13
content:    Support module with sequence utility functions which for some reason
            are missing from Biopython.
'''
# Functions
def expand_ambiguous_seq(seq, seqtype='DNA'):
    '''Expand an ambiguous seq into all possible unambiguous ones'''
    if seqtype == 'DNA':
        from Bio.Data.IUPACData import ambiguous_dna_values as ttable
    elif seqtype == 'RNA':
        from Bio.Data.IUPACData import ambiguous_rna_values as ttable

    # Make a list of possibilities for each character
    amb = (ttable[c] for c in seq)

    # Generate all combinations
    from itertools import product
    seqs = map(''.join, product(*amb))

    return seqs

