# vim: fdm=marker
'''
author:     Fabio Zanini
date:       26/09/13
content:    Manage reference sequences.
'''
# Modules
from Bio import SeqIO

from mapping.filenames import get_HXB2_entire, get_NL43_entire, get_F10_entire, \
        get_HXB2_fragmented, get_NL43_fragmented, get_F10_fragmented



# Functions
def load_HXB2(cropped=False, fragment=None):
    '''Load HXB2 reference sequence'''
    if fragment is None:
        return SeqIO.read(get_HXB2_entire(cropped=cropped), 'fasta')
    else:
        return SeqIO.read(get_HXB2_fragmented(fragment), 'fasta')


def load_NL43(fragment=None):
    '''Load NL4-3 reference sequence'''
    if fragment is None:
        return SeqIO.read(get_NL43_entire(), 'fasta')
    else:
        return SeqIO.read(get_NL43_fragmented(fragment), 'fasta')


def load_F10(fragment=None):
    '''Load F10 reference sequence'''
    if fragment is None:
        return SeqIO.read(get_F10_entire(), 'fasta')
    else:
        return SeqIO.read(get_F10_fragmented(fragment), 'fasta')
