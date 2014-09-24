# vim: fdm=marker
'''
author:     Fabio Zanini
date:       30/01/14
content:    Functions for getting standard filenames.
'''
# Globals
root_data_folder = '/ebio/ag-neher/share/data/PacBio_HIV_Karolinska/'
reference_folder = root_data_folder+'reference/'



# Functions
def get_seqrun_foldername(seq_run):
    '''Get the folder name of a sequencing run'''
    return root_data_folder+seq_run+'/'


def get_custom_reference_filename(reference, format='fasta'):
    '''Get the filename of a custom reference sequence, in one piece'''
    filename = reference
    filename = filename+'.'+format
    return reference_folder+filename


def get_raw_sequence_filename(data_folder, filename):
    '''Get the raw PacBio sequence filename'''
    fn = filename+'_reads.fastq.gz'
    return data_folder+'ccs_reads/'+fn


def get_reference_premap_filename(data_folder, samplename, fragment=None):
    '''Get the filename of the reference used from premapping'''
    fn = 'reference'
    if fragment is not None:
        fn = fn+'_'+fragment
    fn = fn+'.fasta'
    fn = data_folder+samplename+'/premapped/'+fn
    return fn


def get_premapped_filename(data_folder, samplename, type='bam'):
    '''Get the filename of the readed mapped to reference to split into fragments'''
    filename = 'premapped'
    filename = data_folder+samplename+'/premapped/'+filename+'.'+type
    return filename


def get_divided_filenames(data_folder, samplename, fragments, type='bam'):
    '''Get the filenames of the BAM files divided by fragment'''
    filename = 'divided'
    filename = 'divided/'+filename
    filename = data_folder+samplename+'/'+filename
    filenames = []
    for fragment in (list(fragments) + ['ambiguous', 'crossmapped',
                                        'unmapped', 'low_quality']):
        fnf = filename+'_'+fragment+'.'+type
        filenames.append(fnf)
    return filenames


def get_fragment_positions_filename(data_folder, samplename):
    '''Get the filename of the positions of fragments in the reference for premap'''
    filename = 'fragment_positions_premapped.dat'
    filename = 'divided/'+filename
    return data_folder+samplename+'/'+filename




# SUMMARY
def get_divide_summary_filename(data_folder, samplename):
    '''Get the filename of the summary of the division into fragments'''
    filename = 'summary_divide.txt'
    filename = 'divided/'+filename
    filename = data_folder+samplename+'/'+filename
    return filename


