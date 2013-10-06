# vim: fdm=marker
'''
author:     Fabio Zanini
date:       26/09/13
content:    Manage reference sequences.
'''
# Modules
from Bio import SeqIO

from mapping.filenames import get_HXB2_entire



# Functions
def load_HXB2(cropped=False):
    '''Load HXB2 reference sequence'''
    return SeqIO.read(get_HXB2_entire(cropped=cropped), 'fasta')


def load_NL43():
    '''Load NL4-3 reference sequence'''
    return SeqIO.read(get_NL43_entire(), 'fasta')
