# vim: fdm=indent
'''
author:     Fabio Zanini
date:       22/08/13
content:    Module containing all filenames of the analysis in one place.
'''
# Modules
import os

from mapping.adapter_info import foldername_adapter



# Globals
reference_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/reference/'



# Functions
def get_consensus_filename(data_folder, adaID, fragment, subsample=False,
                           trim_primers=False):
    '''Find the filename of the final consensus'''
    filename = 'consensus_'+fragment
    if trim_primers:
        filename = filename+'_trim_primers'
    filename = filename+'.fasta'
    filename = foldername_adapter(adaID)+filename
    if subsample:
        filename = 'subsample/'+filename
    return data_folder+filename


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


def get_mutations_file(data_folder, adaID, fragment, subsample=False):
    '''Get the filename with the mutations for all reads'''
    filename = 'mutations_'+fragment+'.pickle'
    filename = foldername_adapter(adaID)+filename
    if subsample:
        filename = 'subsample/'+filename
    return data_folder+filename


def get_allele_counts_filename(data_folder, adaID, fragment, subsample=False):
    '''Get the filename with the allele counts for all reads'''
    filename = 'allele_counts_'+fragment+'.npy'
    filename = foldername_adapter(adaID)+filename
    if subsample:
        filename = 'subsample/'+filename
    return data_folder+filename


def get_coallele_counts_filename(data_folder, adaID, fragment, subsample=False):
    '''Get the filename with the allele counts for all reads'''
    filename = 'coallele_counts_'+fragment+'.npy'
    filename = foldername_adapter(adaID)+filename
    if subsample:
        filename = 'subsample/'+filename
    return data_folder+filename


def get_insert_counts_filename(data_folder, adaID, fragment, subsample=False):
    '''Get the filename with the insert counts for all reads'''
    filename = 'insert_counts_'+fragment+'.pickle'
    filename = foldername_adapter(adaID)+filename
    if subsample:
        filename = 'subsample/'+filename
    return data_folder+filename


def get_coverage_filename(data_folder, adaID, fragment, subsample=False):
    '''Get the filename with the coverage'''
    filename = 'coverage_'+fragment+'.npy'
    filename = foldername_adapter(adaID)+filename
    if subsample:
        filename = 'subsample/'+filename
    return data_folder+filename


def get_allele_frequencies_filename(data_folder, adaID, fragment, subsample=False):
    '''Get the filename with the corrected allele frequencies'''
    filename = 'allele_frequencies_'+fragment+'.npy'
    filename = foldername_adapter(adaID)+filename
    if subsample:
        filename = 'subsample/'+filename
    return data_folder+filename


def get_HXB2_fragmented(fragment, ext=True):
    '''Get the filename of the reference HXB2 alignment, divided in fragments'''
    filename = 'HXB2_'+fragment
    if ext:
        filename = filename+'.fasta'
    return reference_folder+filename


def get_HXB2_entire(cropped=False):
    '''Get the filename of the reference HXB2 sequence, in one piece'''
    filename = 'HXB2'
    if cropped:
        filename = filename+'_cropped_F1_F6'
    filename = filename+'.fasta'
    return reference_folder+filename


def get_NL43_entire():
    '''Get the filename of the reference NL4-3 sequence, in one piece'''
    filename = 'NL4-3'
    filename = filename+'.fasta'
    return reference_folder+filename


def get_HXB2_index_file(fragment='F0', ext=True):
    '''Get the index filename, with or w/o extension'''
    filename = 'HXB2'
    if fragment != 'F0':
        filename = filename+'_'+fragment
    filename = reference_folder+'hash/'+filename
    if ext:
        filename = filename+'.stidx'
    return filename


def get_HXB2_hash_file(fragment='F0', ext=True):
    '''Get the index filename, with or w/o extension'''
    filename = 'HXB2'
    if fragment != 'F0':
        filename = filename+'_'+fragment
    filename = reference_folder+'hash/'+filename
    if ext:
        filename = filename+'.sthash'
    return filename


def get_read_filenames(data_folder, adaID, fragment=None, subsample=False,
                       filtered=True, premapped=False, suffix=''):
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
                fn = [fn+'_F'+str(j)+suffix+'.fastq' for j in xrange(1, 7)] +\
                     [fn+'_ambiguous.fastq'] +\
                     [fn+'_unmapped.fastq']
            else:
                fn = fn+'_'+fragment+suffix+'.fastq'
        else:
            fn = fn+suffix+'.fastq'
        filenames[i] = fn
    return filenames


def get_premapped_file(data_folder, adaID, type='bam', subsample=False, bwa=False):
    '''Get the filename of the readed mapped to HXB2 to split into fragments'''
    filename = 'premapped_to_HXB2'
    if bwa:
        filename = filename + '_bwa'

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


def get_mapped_filename(data_folder, adaID, fragment, type='bam', subsample=False,
                        bwa=False, filtered=False, sort=False):
    '''Get the filename of the mapped reads onto consensus'''
    filename = fragment
    if bwa:
        filename = filename + '_bwa'
    if filtered:
        filename = filename + '_filtered'
    if sort:
        filename = filename + '_sorted'

    if type == 'sam':
        filename = filename + '.sam'
    elif type == 'bam':
        filename = filename + '.bam'
    else:
        raise ValueError('Type of mapped reads file not recognized')
    filename = 'mapped/'+filename
    filename = foldername_adapter(adaID)+filename
    if subsample:
        filename = 'subsample/'+filename
    return data_folder+filename


def get_raw_read_files(dataset):
    '''Get the raw files which we obtain from Joerg/Xi'''
    data_folder = dataset['folder'].rstrip('/')+'/'
    return {key: data_folder+subdir
            for key, subdir in dataset['raw_data'].iteritems()}


def get_read_unpaired_filename(data_folder, adaID, subsample=False):
    '''Get the reads pairs for which one read is low quality'''
    fn = 'reads_unpaired.fastq'
    fn = foldername_adapter(adaID)+fn
    if subsample:
        fn = 'subsample/'+fn
    fn = data_folder+fn
    return fn


def get_phix_filename():
    '''Get the phiX sequence filename'''
    filename = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/reference/phiX_genome.fasta'
    return filename


def get_mapped_phix_filename(data_folder, type='bam', filtered=False, sort=False):
    '''Get the filename of the mapped reads onto PhiX'''
    filename = 'mapped_to_phix'
    if filtered:
        filename = filename + '_filtered'
    if sort:
        filename = filename + '_sorted'
    if type == 'sam':
        filename = filename + '.sam'
    elif type == 'bam':
        filename = filename + '.bam'
    else:
        raise ValueError('Type of mapped reads file not recognized')
    return data_folder+'phix/'+filename


def get_allele_counts_phix_filename(data_folder):
    '''Get the filename of the allele counts of PhiX'''
    return data_folder+'phix/'+'allele_counts.npy'


def get_insert_counts_phix_filename(data_folder):
    '''Get the filename of the insert counts of PhiX'''
    return data_folder+'phix/'+'insert_counts.pickle'


def get_unclassified_reads_filenames(data_folder, filtered=False):
    '''Get the filenames of the unclassified reads'''
    filenames = ['read1', 'read2', 'adapter']
    if filtered:
        filenames = [f+'_filtered_trimmed' if 'read' in f else f for f in filenames]
    filenames = [data_folder+'unclassified_reads/'+f+'.fastq' for f in filenames]
    return filenames
