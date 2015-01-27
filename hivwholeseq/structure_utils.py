# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/01/15
content:    Support module for structural computations.
'''
# Functions
def get_resname1(res, throw=False):
    '''Get residue name in one letter format'''
    from Bio.Data.IUPACData import protein_letters_3to1_extended as table
    rname = res.get_resname().capitalize()
    if rname not in table:
        if throw:
            raise ValueError('Residue '+rname+' not found in conversion table')
        else:
            return 'X'

    return table[rname]


def get_chainseq(chain, **kwargs):
    '''Get the sequence of a PDB chain'''
    return ''.join(get_resname1(res, **kwargs) for res in chain.get_residues())
