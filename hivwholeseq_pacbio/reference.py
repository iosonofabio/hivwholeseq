# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/09/14
content:    Support module for reference sequences.
'''
# Modules
from Bio import SeqIO
from hivwholeseq_pacbio.sequencing.filenames import get_custom_reference_filename
from hivwholeseq.sequence_utils import correct_genbank_features_load, \
        correct_genbank_features_save



# Functions
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

