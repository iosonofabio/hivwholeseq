# vim: fdm=indent
'''
author:     Fabio Zanini
date:       22/08/13
content:    Module containing all filenames of the analysis in one place.
'''
# Modules
import os

from mapping.adapter_info import foldername_adapter


# Functions
def get_last_consensus_number(data_folder, adaID):
    '''Find the number of the final consensus after iterative mapping'''
    # Directory to read
    dirname = foldername_adapter(adaID)

    # Find the last (fragmented) consensus
    g = os.walk(data_folder+'subsample/'+dirname)
    fns = g.next()[2]
    fns = filter(lambda x: ('consensus' in x) and ('fragmented.fasta' in x), fns)
    consensus_numbers = map(lambda x: int(x.split('_')[1]), fns)
    cons_max = max(consensus_numbers)
    return cons_max


def get_last_reference(data_folder, adaID, ext=True):
    '''Find the filename of the last consensus after iterative mapping'''
    cons_max = get_last_consensus_number(data_folder, adaID)
    filename = 'consensus_'+str(cons_max)+'_fragmented'

    if ext:
        filename = filename + '.fasta'
    return data_folder+'subsample/'+foldername_adapter(adaID)+filename


def get_last_mapped(data_folder, adaID, type='bam'):
    '''Find the filename of the mapped reads to the last consensus'''
    cons_max = get_last_consensus_number(data_folder, adaID)
    filename = 'mapped_to_consensus_'+str(cons_max)+'_fragmented'
    if type == 'sam':
        filename = filename + '.sam'
    elif type == 'bam':
        filename = filename + '.bam'
    else:
        raise ValueError('Type of mapped reads file not recognized')
    return data_folder+foldername_adapter(adaID)+filename


def get_mutations_file(data_folder, adaID):
    '''Get the filename with the mutations for all reads'''
    filename = 'mutations.pickle'
    filename = data_folder+foldername_adapter(adaID)+filename
    return filename
