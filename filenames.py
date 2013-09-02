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
def get_consensus_filename(data_folder, adaID, fragment):
    '''Find the filename of the final consensus'''
    filename = 'consensus_'+fragment+'.fasta'
    filename = data_folder+foldername_adapter(adaID)+filename
    return filename


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


def get_last_mapped(data_folder, adaID, type='bam', filtered=False):
    '''Find the filename of the mapped reads to the last consensus'''
    cons_max = get_last_consensus_number(data_folder, adaID)
    filename = 'mapped_to_consensus_'+str(cons_max)+'_fragmented'
    if filtered:
        filename = filename + '_filtered'
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
    return data_folder+foldername_adapter(adaID)+filename


def get_allele_counts_filename(data_folder, adaID):
    '''Get the filename with the allele counts for all reads'''
    filename = 'allele_counts.npy'
    return data_folder+foldername_adapter(adaID)+filename


def get_insert_counts_filename(data_folder, adaID):
    '''Get the filename with the insert counts for all reads'''
    filename = 'insert_counts.pickle'
    return data_folder+foldername_adapter(adaID)+filename


def get_coverage_filename(data_folder, adaID):
    '''Get the filename with the coverage'''
    filename = 'coverage.npy'
    return data_folder+foldername_adapter(adaID)+filename


def get_HXB2_fragmented(data_folder, fragment):
    '''Get the filename of the reference HXB2 alignment, divided in fragments'''
    filename = 'HXB2_'+fragment+'.fasta'
    return data_folder+'reference/'+filename


def get_HXB2_entire(data_folder, cropped=False):
    '''Get the filename of the reference HXB2 alignment, in one piece'''
    filename = 'HXB2'
    if cropped:
        filename = filename+'_cropped_F1_F6'
    filename = filename+'.fasta'
    return data_folder+'reference/'+filename


def get_HXB2_index_file(data_folder, fragment='F0', ext=True):
    '''Get the index filename, with or w/o extension'''
    filename = 'HXB2'
    if fragment != 'F0':
        filename = filename+'_'+fragment
    filename = data_folder+'reference/'+filename
    if ext:
        filename = filename+'.stidx'
    return filename


def get_HXB2_hash_file(data_folder, fragment='F0', ext=True):
    '''Get the index filename, with or w/o extension'''
    filename = 'HXB2'
    if fragment != 'F0':
        filename = filename+'_'+fragment
    filename = data_folder+'reference/'+filename
    if ext:
        filename = filename+'.sthash'
    return filename


def get_read_filenames(data_folder, adaID, subsample=False,
                       filtered=True, premapped=False, fragment=None):
    '''Get the filenames of the demultiplexed reads'''
    filenames = ['read1', 'read2']
    for i,fn in enumerate(filenames):
        if premapped:
            fn = 'premapped/'+fn
        elif filtered:
            fn = fn+'_filtered_trimmed'
        fn = foldername_adapter(adaID)+fn
        if subsample:
            fn = 'subsample/'+fn
        fn = data_folder+fn
        
        # If there was a premap, 6 files have to be made for each orientation
        if premapped:
            if fragment is None:
                fn = [fn+'_F'+str(j)+'.fastq' for j in xrange(1, 7)] + [fn+'_unmapped.fastq']
            else:
                fn = fn+'_'+fragment+'.fastq'
        else:
            fn = fn+'.fastq'
        filenames[i] = fn
    return filenames


def get_premapped_file(data_folder, adaID, type='bam', subsample=False):
    '''Get the filename of the readed mapped to HXB2 to split into fragments'''
    filename = 'premapped_to_HXB2'
    if type == 'sam':
        filename = filename + '.sam'
    elif type == 'bam':
        filename = filename + '.bam'
    else:
        raise ValueError('Type of mapped reads file not recognized')
    filename = 'premapped/'+filename
    filename = foldername_adapter(adaID)+filename
    if subsample:
        filename = 'subsample/'+filename
    return data_folder+filename

