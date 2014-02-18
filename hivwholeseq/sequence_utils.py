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


def pretty_print_pairwise_ali(ali, name1='', name2=''):
    '''Pretty print function for pairwise alignments'''
    from itertools import izip

    for i in xrange(len(ali[0]) / 51 + 1):
        ali1_t = str(ali[0, i * 51: (i+1) * 51].seq)
        ali2_t = str(ali[1, i * 51: (i+1) * 51].seq)
        match_t = []
        for (a1, a2) in izip(ali1_t, ali2_t):
            if a1 == a2:
                match_t.append(' ')
            else:
                match_t.append('x')
        match_t = ''.join(match_t)

        lh = min(max(map(len, [name1, name2])), 6)
        print name1[:lh]+(' ' * max(0, lh - len(name1)))+':', ali1_t
        print (' ' * (lh + 1)), match_t
        print name2[:lh]+(' ' * max(0, lh - len(name2)))+':', ali2_t
        print

