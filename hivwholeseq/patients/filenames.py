# vim: fdm=marker
'''
author:     Fabio Zanini
date:       25/10/13
content:    Module containing all filenames for the patient analysis in one place.
'''
# Modules



# Globals
# FIXME
root_data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/'



# Functions
def get_foldername(pname, root_data_folder=root_data_folder):
    '''Get the folder name of the data from a patient'''
    foldername = 'patients/'+pname+'/'
    foldername = root_data_folder+foldername
    return foldername


def get_initial_consensus_filename(pname, fragment,
                                   root_data_folder=root_data_folder,
                                   format='fasta'):
    '''Get the filename of the initial consensus for a patient'''
    filename = 'consensus_initial_'+fragment+'.'+format
    filename = get_foldername(pname, root_data_folder=root_data_folder)+filename
    return filename


def get_mapped_to_initial_filename(pname, samplename, fragment, type='bam',
                                   part=None, unsorted=False,
                                   filtered=False,
                                   root_data_folder=root_data_folder):
    '''Get the filename of the mapped reads to initial consensus'''
    filename = fragment
    if part is not None:
        filename = filename+'_part'+str(part)
    elif unsorted:
        filename = filename+'_unsorted'
    if filtered:
        filename = filename+'_filtered'
    filename = filename+'.'+type
    filename = samplename+'/mapped_to_initial/'+filename
    filename = get_foldername(pname, root_data_folder=root_data_folder)+filename
    return filename


def get_initial_index_filename(pname, fragment, ext=True,
                               root_data_folder=root_data_folder):
    '''Get the index filename, with or w/o extension'''
    filename = 'consensus_initial_'+fragment
    filename = 'hash/'+filename
    filename = get_foldername(pname, root_data_folder=root_data_folder)+filename
    if ext:
        filename = filename+'.stidx'
    return filename


def get_initial_hash_filename(pname, fragment, ext=True,
                              root_data_folder=root_data_folder):
    '''Get the index filename, with or w/o extension'''
    filename = 'consensus_initial_'+fragment
    filename = 'hash/'+filename
    filename = get_foldername(pname, root_data_folder=root_data_folder)+filename
    if ext:
        filename = filename+'.sthash'
    return filename


def get_consensi_alignment_filename(pname, fragment,
                                    root_data_folder=root_data_folder):
    '''Get the MSA of all consensi of the patient, sorted by time point'''
    filename = 'consensi_alignment_'+fragment+'.fasta'
    filename = get_foldername(pname, root_data_folder=root_data_folder)+filename
    return filename


def get_consensi_alignment_genomewide_filename(pname,
                                               root_data_folder=root_data_folder):
    '''Get the MSA of all consensi of the patient, sorted by time point'''
    filename = 'consensi_alignment_genomewide.fasta'
    filename = get_foldername(pname, root_data_folder=root_data_folder)+filename
    return filename


def get_allele_count_trajectories_filename(pname,
                                           fragment,
                                           root_data_folder=root_data_folder):
    '''Get the matrix with allele counts on the initial consensus'''
    filename = 'allele_frequency_counts'+fragment+'.npy'
    filename = get_foldername(pname, root_data_folder=root_data_folder)+filename
    return filename


def get_allele_frequency_trajectories_filename(pname,
                                               fragment,
                                               root_data_folder=root_data_folder):
    '''Get the matrix with allele frequencies on the initial consensus'''
    filename = 'allele_frequency_trajectories_'+fragment+'.npy'
    filename = get_foldername(pname, root_data_folder=root_data_folder)+filename
    return filename


def get_allele_cocounts_filename(pname, samplename, fragment,
                                 root_data_folder=root_data_folder):
    '''Get the matrix of allele cocounts on the initial consensus'''
    filename = 'allele_cocounts_'+fragment+'.npy'
    filename = get_foldername(pname, root_data_folder=root_data_folder)+samplename+'/'+filename
    return filename


# FIGURES
def get_figure_folder(pname):
    '''Get the folder for figures for this sample'''
    folder = get_foldername(pname, root_data_folder=root_data_folder)+'figures/'
    return folder


def get_allele_frequency_trajectory_figure_filename(pname, fragment_or_gene, format='png'):
    '''Get the filename of the plot of allele frequency trajectories'''
    folder = get_figure_folder(pname)
    fn = folder+'allele_freq_traj_'+fragment_or_gene+'.'+format
    return fn


# SUMMARY
def get_map_initial_summary_filename(pname, samplename, fragment,
                                     root_data_folder=root_data_folder):
    '''Get the filename of the summary of the division into fragments'''
    filename = 'summary_map_initial_'+fragment+'.txt'
    filename = samplename+'/mapped_to_initial/'+filename
    filename = get_foldername(pname, root_data_folder=root_data_folder)+filename
    return filename


def get_filter_mapped_init_summary_filename(pname, samplename, fragment,
                                     root_data_folder=root_data_folder):
    '''Get the filename of the summary of the post-map filtering'''
    filename= 'summary_filter_mapped_init_'+fragment+'.txt'
    filename = samplename+'/mapped_to_initial/'+filename
    filename = get_foldername(pname, root_data_folder=root_data_folder)+filename
    return filename
