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
def load_HXB2(data_folder, cropped=False):
    '''Load HXB2 reference sequence'''
    return SeqIO.read(get_HXB2_entire(data_folder, cropped=cropped), 'fasta')


def load_NL43(data_folder):
    '''Load NL4-3 reference sequence'''
    return SeqIO.read(get_NL43_entire(data_folder, cropped=cropped), 'fasta')
