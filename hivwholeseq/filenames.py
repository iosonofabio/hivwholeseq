# vim: fdm=indent
'''
author:     Fabio Zanini
date:       22/08/13
content:    Module containing all filenames of the analysis in one place.
'''
# Modules
import os

from hivwholeseq.adapter_info import foldername_adapter
import hivwholeseq



# Globals
root_data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/'
reference_folder = root_data_folder+'reference/'
table_filename = hivwholeseq.__path__[0] + '/tables/HIV_reservoir_all.xlsx'



# Functions
def get_seqrun_foldername(seq_run):
    '''Get the path of the folder for a sequencing run'''
    fn = root_data_folder+seq_run+'/'
    return fn


def get_reference_premap_filename(data_folder, adaID, fragment=None):
    '''Get the filename of the reference used from premapping'''
    fn = 'reference'
    if fragment is not None:
        fn = fn+'_'+fragment
    fn = fn+'.fasta'
    fn = data_folder+foldername_adapter(adaID)+'premapped/'+fn
    return fn


def get_reference_premap_index_filename(data_folder, adaID, ext=True):
    '''Get the filename of the stampy index of the reference used for premapping'''
    fn = 'reference'
    if ext:
        fn = fn + '.stidx'
    fn = data_folder+foldername_adapter(adaID)+'premapped/'+fn
    return fn


def get_reference_premap_hash_filename(data_folder, adaID, ext=True):
    '''Get the filename of the stampy hash of the reference used for premapping'''
    fn = 'reference'
    if ext:
        fn = fn + '.sthash'
    fn = data_folder+foldername_adapter(adaID)+'premapped/'+fn
    return fn


def get_consensus_filename(data_folder, adaID=None, fragment=None, trim_primers=True):
    '''Find the filename of the final consensus'''
    filename = 'consensus'
    if fragment is not None:
        filename = filename+'_'+fragment
    if not trim_primers:
        filename = filename+'_with_primers'
    filename = filename+'.fasta'
    if adaID is not None:
        filename = foldername_adapter(adaID)+filename
    return data_folder+filename


def get_consensus_old_filename(data_folder, adaID, fragment, trim_primers=True):
    '''Find the filename of the final consensus'''
    filename = 'consensus_old_'+fragment
    if not trim_primers:
        filename = filename+'_with_primers'
    filename = filename+'.fasta'
    filename = foldername_adapter(adaID)+filename
    return data_folder+filename


def get_mutations_file(data_folder, adaID, fragment):
    '''Get the filename with the mutations for all reads'''
    filename = 'mutations_'+fragment+'.pickle'
    filename = foldername_adapter(adaID)+filename
    return data_folder+filename


def get_allele_counts_filename(data_folder, adaID, fragment):
    '''Get the filename with the allele counts for all reads'''
    filename = 'allele_counts_'+fragment+'.npy'
    filename = foldername_adapter(adaID)+filename
    return data_folder+filename


def get_coallele_counts_filename(data_folder, adaID, fragment):
    '''Get the filename with the allele counts for all reads'''
    filename = 'coallele_counts_'+fragment+'.npy'
    filename = foldername_adapter(adaID)+filename
    return data_folder+filename


def get_insert_counts_filename(data_folder, adaID, fragment):
    '''Get the filename with the insert counts for all reads'''
    filename = 'insert_counts_'+fragment+'.pickle'
    filename = foldername_adapter(adaID)+filename
    return data_folder+filename


def get_coverage_filename(data_folder, adaID, fragment):
    '''Get the filename with the coverage'''
    filename = 'coverage_'+fragment+'.npy'
    filename = foldername_adapter(adaID)+filename
    return data_folder+filename


def get_allele_frequencies_filename(data_folder, adaID, fragment):
    '''Get the filename with the corrected allele frequencies'''
    filename = 'allele_frequencies_'+fragment+'.npy'
    filename = foldername_adapter(adaID)+filename
    return data_folder+filename


def get_HXB2_fragmented(fragment, ext=True, trim_primers=False):
    '''Get the filename of the reference HXB2 sequence, divided in fragments'''
    filename = 'HXB2_'+fragment
    if trim_primers:
        filename = filename+'_trim_primers'
    if ext:
        filename = filename+'.fasta'
    return reference_folder+filename


def get_NL43_fragmented(fragment, ext=True, trim_primers=False):
    '''Get the filename of the reference NL4-3 sequence, divided in fragments'''
    filename = 'NL4-3_'+fragment
    if trim_primers:
        filename = filename+'_trim_primers'
    if ext:
        filename = filename+'.fasta'
    return reference_folder+filename


def get_F10_fragmented(fragment, ext=True, trim_primers=False):
    '''Get the filename of the reference F10 sequence, divided in fragments'''
    filename = 'F10_'+fragment
    if trim_primers:
        filename = filename+'_trim_primers'
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


def get_F10_entire():
    '''Get the filename of the reference F10 (ZM246F) sequence, in one piece'''
    filename = 'F10'
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


def get_read_filenames(data_folder, adaID=None, fragment=None, suffix='',
                       gzip=False, trimmed=False):
    '''Get the filenames of the demultiplexed reads'''
    filenames = ['read1', 'read2']
    for i,fn in enumerate(filenames):
        if adaID is not None:
            fn = foldername_adapter(adaID)+fn
        fn = data_folder+fn
        if trimmed:
            fn = fn+'_trimmed'
        fn = fn+suffix+'.fastq'
        if gzip:
            fn = fn+'.gz' 
        filenames[i] = fn
    return filenames


def get_premapped_filename(data_folder, adaID=None, type='bam', bwa=False,
                           part=None, unsorted=False):
    '''Get the filename of the readed mapped to reference to split into fragments'''
    filename = 'premapped'
    filename = 'premapped/'+filename
    if adaID is not None:
        filename = foldername_adapter(adaID)+filename
    if part is not None:
        filename = filename+'_part'+str(part)
    elif unsorted:
        filename = filename+'_unsorted'

    if bwa:
        filename = filename + '_bwa'
    if type == 'sam':
        filename = filename + '.sam'
    elif type == 'bam':
        filename = filename + '.bam'
    else:
        raise ValueError('Type of mapped reads file not recognized')

    return data_folder+filename


def get_fragment_positions_filename(data_folder, adaID):
    '''Get the filename of the positions of fragments in the reference for premap'''
    filename = 'fragment_positions_premapped.dat'
    filename = 'divided/'+filename
    return data_folder+foldername_adapter(adaID)+filename


def get_divided_filenames(data_folder, adaID=None, fragments=None, type='bam'):
    '''Get the filenames of the BAM files divided by fragment'''
    filename = 'divided'
    filename = 'divided/'+filename
    if adaID is not None:
        filename = foldername_adapter(adaID)+filename
    filename = data_folder+filename
    filenames = []
    for fragment in (list(fragments) + ['ambiguous', 'crossmapped',
                                        'unmapped', 'low_quality']):
        fnf = filename+'_'+fragment+'.'+type
        filenames.append(fnf)
    return filenames


def get_divided_filename(data_folder, adaID=None, fragment=None, type='bam', chunk=None):
    '''Get the filename of the BAM files divided for a single fragment'''
    filename = 'divided'
    filename = 'divided/'+filename
    if adaID is not None:
        filename = foldername_adapter(adaID)+filename
    filename = data_folder+filename
    filename = filename+'_'+fragment
    if chunk is not None:
        filename = filename+'_chunk_'+str(chunk)
    filename = filename+'.'+type
    return filename


def get_mapped_filename(data_folder, adaID=None, fragment=None, type='bam', 
                        bwa=False, filtered=False, sort=False, part=None, unsorted=False,
                        rescue=False):
    '''Get the filename of the mapped reads onto consensus'''
    if fragment is None:
        raise ValueError('Select a fragment')
    filename = fragment
    if rescue:
        filename = filename + '_rescue'
    if bwa:
        filename = filename + '_bwa'
    if filtered:
        filename = filename + '_filtered'
    if sort:
        filename = filename + '_sorted'
    elif part is not None:
        filename = filename+'_part'+str(part)
    elif unsorted:
        filename = filename+'_unsorted'

    filename = 'mapped/'+filename+'.'+type
    if adaID is not None:
        filename = foldername_adapter(adaID)+filename
    return data_folder+filename


def get_mapped_suspicious_filename(data_folder, adaID, fragment, type='bam'):
    '''The the filename of the mapped reads with many mutations from consensus'''
    filename = fragment+'_suspicious.'+type
    filename = data_folder+foldername_adapter(adaID)+'mapped/'+filename
    return filename


def get_raw_read_files(dataset):
    '''Get the raw files which we obtain from Joerg/Xi'''
    data_folder = dataset['folder'].rstrip('/')+'/'
    return {key: data_folder+subdir
            for key, subdir in dataset['raw_data'].iteritems()}


def get_read_unpaired_filename(data_folder, adaID):
    '''Get the reads pairs for which one read is low quality'''
    fn = 'reads_unpaired.fastq'
    fn = foldername_adapter(adaID)+fn
    fn = data_folder+fn
    return fn


def get_phix_filename():
    '''Get the phiX sequence filename'''
    filename = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/reference/phiX_genome.fasta'
    return filename


def get_mapped_phix_filename(data_folder, type='bam', filtered=False, sort=False,
                             part=None, unsorted=False):
    '''Get the filename of the mapped reads onto PhiX'''
    filename = 'mapped_to_phix'

    if part is not None:
        filename = filename+'_part'+str(part)
    elif unsorted:
        filename = filename+'_unsorted'

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


def get_consensus_phix_filename(data_folder):
    '''Get the filename of the consensus of PhiX'''
    return data_folder+'phix/'+'consensus.fasta'


def get_allele_counts_phix_filename(data_folder, qual_min=None):
    '''Get the filename of the allele counts of PhiX'''
    fn = data_folder+'phix/'+'allele_counts'
    if qual_min is not None:
        fn = fn+'_qualmin_'+str(qual_min)
    fn = fn+'.npy'
    return fn


def get_insert_counts_phix_filename(data_folder, qual_min=None):
    '''Get the filename of the insert counts of PhiX'''
    fn = data_folder+'phix/'+'insert_counts'
    if qual_min is not None:
        fn = fn+'_qualmin_'+str(qual_min)
    fn = fn+'.pickle'
    return fn


def get_unclassified_reads_filenames(data_folder, filtered=False, gzip=False):
    '''Get the filenames of the unclassified reads'''
    filenames = ['read1', 'read2', 'adapter']
    if filtered:
        filenames = [f+'_filtered_trimmed' if 'read' in f else f for f in filenames]
    filenames = [data_folder+'unclassified_reads/'+f+'.fastq' for f in filenames]
    if gzip:
        filenames = [f+'.gz' for f in filenames]
    return filenames


def get_merged_consensus_filename(data_folder, adaID=None,
                                  fragments=['F1', 'F2', 'F3', 'F4', 'F5', 'F6']):
    '''Get the merged consensus of several fragments'''
    filename = 'consensus_'+'-'.join(fragments)+'.fasta'
    if adaID is not None:
        filename = foldername_adapter(adaID)+filename
    filename = data_folder+filename
    return filename


def get_merged_allele_frequencies_filename(data_folder, adaID,
                                    fragments=['F1', 'F2', 'F3', 'F4', 'F5', 'F6']):
    '''Get the merged allele frequencies of several fragments'''
    filename = 'allele_frequencies_'+'-'.join(fragments)+'.fasta'
    filename = data_folder+foldername_adapter(adaID)+filename
    return filename


def get_custom_reference_filename(reference, format='fasta'):
    '''Get the filename of a custom reference sequence, in one piece'''
    filename = reference
    filename = filename+'.'+format
    return reference_folder+filename


def get_custom_index_filename_fun(reference):
    '''Get the stampy index filename of a custom reference sequence'''
    filename = reference
    filename = reference_folder+'hash/'+filename
    fun = lambda ext: [filename+'.stidx' if ext else filename][0]
    return fun


def get_custom_hash_filename_fun(reference):
    '''Get the stampy index filename of a custom reference sequence'''
    filename = reference
    filename = reference_folder+'hash/'+filename
    fun = lambda ext=True: [filename+'.sthash' if ext else filename][0]
    return fun


# FIGURES
def get_figure_folder(data_folder, adaID=None):
    '''Get the folder for figures for this sample'''
    folder = 'figures/'
    if adaID is not None:
        folder = foldername_adapter(adaID)+folder
    folder = data_folder+folder
    return folder


def get_quality_along_reads_filename(data_folder, adaID=None, ext='png', simple=False):
    '''Get the filename of the quality along the reads'''
    filename = 'quality_along_reads'
    if simple:
        filename = filename+'_simple'
    filename = filename+'.'+ext
    filename = get_figure_folder(data_folder, adaID)+filename
    return filename


def get_insert_size_distribution_cumulative_filename(data_folder, adaID, fragment,
                                                     ext='png'):
    '''Get filename of the cumulative distribution of insert sizes'''
    filename = 'insert_size_cumulative_distribution'
    if fragment == 'premapped':
        filename = filename+'_premapped'
    else:
        filename = filename+'_mapped_'+fragment
    filename = filename+'.'+ext
    filename = get_figure_folder(data_folder, adaID)+filename
    return filename


def get_insert_size_distribution_filename(data_folder, adaID, fragment,
                                          ext='png'):
    '''Get filename of the distribution of insert sizes'''
    filename = 'insert_size_distribution'
    if fragment == 'premapped':
        filename = filename+'_premapped'
    else:
        filename = filename+'_mapped_'+fragment
    filename = filename+'.'+ext
    filename = get_figure_folder(data_folder, adaID)+filename
    return filename


def get_coverage_figure_filename(data_folder, adaID, fragment, ext='png'):
    '''Get the filename of the coverage report figure'''
    filename = 'coverage'
    if fragment == 'premapped':
        filename = filename+'_premapped'
    else:
        filename = filename+'_mapped_'+fragment
    filename = filename+'.'+ext
    filename = get_figure_folder(data_folder, adaID)+filename
    return filename
    

def get_overlap_nu_figure_filename(data_folder, adaID, fragments, ext='png'):
    '''Get the filename of the coverage report figure'''
    filename = 'overlap_nu_'+fragments
    filename = 'overlap/'+filename
    filename = filename+'.'+ext
    filename = get_figure_folder(data_folder, adaID)+filename
    return filename


def get_distance_from_consensus_figure_filename(data_folder, adaID, fragment,
                                                yscale='linear',
                                                sliding_window=False,
                                                cumulative=False,
                                                ext='png'):
    '''Get the filename of the figure of distance from fragment consensus'''
    filename = 'distance_from_consensus_histogram_'+fragment
    if sliding_window:
        filename = filename+'_sliding_window'
    if cumulative:
        filename = filename+'_cumulative'
    if yscale != 'linear':
        filename = filename+'_y'+yscale
    filename = filename+'.'+ext
    filename = get_figure_folder(data_folder, adaID)+filename
    return filename


def get_SFS_figure_filename(data_folder, adaID, fragment,
                            yscale='linear',
                            cumulative=False,
                            ext='png'):
    '''Get the filename of the SFS figure'''
    filename = 'SFS_'+fragment
    if cumulative:
        filename = filename+'_cumulative'
    if yscale != 'linear':
        filename = filename+'_y'+yscale
    filename = filename+'.'+ext
    filename = get_figure_folder(data_folder, adaID)+filename
    return filename


def get_minor_allele_frequency_figure_filename(data_folder, adaID, fragments,
                                               only_filtered=False,
                                               ext='png'):
    '''Get the filename of the figure of the minor allele frequency along the genome'''
    filename = 'minor_nu_'+'_'.join(fragments)
    if only_filtered:
        filename = filename+'_filtered'
    filename = filename+'.'+ext
    filename = get_figure_folder(data_folder, adaID)+filename
    return filename
    

def get_minor_allele_frequency_merged_figure_filename(data_folder, adaID, fragments,
                                               ext='png'):
    '''Get the filename of the figure of the minor allele frequency along the genome'''
    filename = 'minor_nu_merged'+'_'.join(fragments)
    filename = filename+'_filtered'
    filename = filename+'.'+ext
    filename = get_figure_folder(data_folder, adaID)+filename
    return filename


def get_consensi_alignment_dataset_filename(data_folder, fragment):
    '''Get the filename of the alignment of consensi from a dataset'''
    fn = 'ali_consensi_'+fragment+'.fasta'
    fn = data_folder+fn
    return fn

    

# SUMMARY
def get_demultiplex_summary_filename(data_folder):
    '''Get the filename of the summary of demultiplex'''
    filename = 'summary_demultiplex.txt'
    filename = data_folder+filename
    return filename


def get_trim_summary_filename(data_folder, adaID):
    '''Get the filename of the trim low quality'''
    filename = 'summary_trim.txt'
    filename = data_folder+foldername_adapter(adaID)+filename
    return filename


def get_premap_summary_filename(data_folder, adaID):
    '''Get the filename of the premap to reference'''
    filename = 'summary_premapped.txt'
    filename = 'premapped/'+filename
    filename = data_folder+foldername_adapter(adaID)+filename
    return filename


def get_divide_summary_filename(data_folder, adaID):
    '''Get the filename of the summary of the division into fragments'''
    filename = 'summary_divide.txt'
    filename = 'divided/'+filename
    filename = data_folder+foldername_adapter(adaID)+filename
    return filename


def get_build_consensus_summary_filename(data_folder, adaID, fragment='general',
                                         iterative=True):
    '''Get the filename of the summary of the iterative consensus'''
    filename = 'summary_build_consensus_'+fragment+'.txt'
    if iterative:
        filename = 'map_iter/'+filename
    filename = data_folder+foldername_adapter(adaID)+filename
    return filename


def get_map_summary_filename(data_folder, adaID, fragment, rescue=False):
    '''Get the filename of the summary of the division into fragments'''
    filename = 'summary_map'+fragment
    if rescue:
        filename = filename+'_rescue'
    filename = filename+'.txt'
    filename = 'mapped/'+filename
    filename = data_folder+foldername_adapter(adaID)+filename
    return filename


def get_filter_mapped_summary_filename(data_folder, adaID, fragment):
    '''Get the filename of the summary of the division into fragments'''
    filename = 'summary_filter_'+fragment+'.txt'
    filename = 'mapped/'+filename
    filename = data_folder+foldername_adapter(adaID)+filename
    return filename


