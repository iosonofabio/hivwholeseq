# vim: fdm=marker
'''
author:     Fabio Zanini
date:       26/09/13
content:    Manage reference sequences.
'''
# Modules
from Bio import SeqIO

from hivwholeseq.filenames import get_HXB2_entire, get_NL43_entire, get_F10_entire, \
        get_HXB2_fragmented, get_NL43_fragmented, get_F10_fragmented, \
        get_custom_reference_filename
from hivwholeseq.sequence_utils import correct_genbank_features_load, \
        correct_genbank_features_save



# Functions
def load_HXB2(cropped=False, fragment=None, trim_primers=False):
    '''Load HXB2 reference sequence'''
    if fragment is None:
        return SeqIO.read(get_HXB2_entire(cropped=cropped), 'fasta')
    else:
        return SeqIO.read(get_HXB2_fragmented(fragment,
                                              trim_primers=trim_primers),
                          'fasta')


def load_NL43(fragment=None, trim_primers=False):
    '''Load NL4-3 reference sequence'''
    if fragment is None:
        return SeqIO.read(get_NL43_entire(), 'fasta')
    else:
        return SeqIO.read(get_NL43_fragmented(fragment,
                                              trim_primers=trim_primers),
                          'fasta')


def load_F10(fragment=None):
    '''Load F10 reference sequence'''
    if fragment is None:
        return SeqIO.read(get_F10_entire(), 'fasta')
    else:
        return SeqIO.read(get_F10_fragmented(fragment,
                                             trim_primers=trim_primers),
                          'fasta')


def load_custom_reference(reference, format='fasta'):
    '''Load a custom reference'''
    record = SeqIO.read(get_custom_reference_filename(reference, format=format), format)

    # BUG: feature id is lost during write, fake it with the 'note' qualifier
    if format in ['gb' , 'genbank']:
        correct_genbank_features_load(record)
    return record


def save_custom_reference(record, reference, format='fasta', molecule='DNA'):
    '''Save a custom reference'''
    if format in ['gb' , 'genbank']:
        correct_genbank_features_save(record, molecule='DNA')

    return SeqIO.write(record,
                       get_custom_reference_filename(reference, format=format),
                       format)

