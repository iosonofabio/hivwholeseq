# vim: fdm=marker
'''
author:     Fabio Zanini
date:       25/10/13
content:    Module containing all filenames for the patient analysis in one place.
'''
# Modules
from hivwholeseq.filenames import root_data_folder


# Globals
root_patient_folder = root_data_folder+'patients/'



# Functions
def get_foldername(pname):
    '''Get the folder name of the data from a patient'''
    foldername = root_patient_folder+pname+'/'
    return foldername


def get_sample_foldername(pname, samplename, PCR=None):
    '''Get the folder name of the data from a patient sample'''
    fn = get_foldername(pname)+'samples/'+samplename+'/'
    if PCR is not None:
        fn = fn+'PCR'+str(PCR)+'/'
    return fn


def get_initial_reference_foldername(pname):
    '''Get the folder name of the initial references'''
    fn = get_foldername(pname)+'reference/'
    return fn


def get_initial_reference_filename(pname, fragment, format='fasta'):
    '''Get the filename of the initial reference for a patient'''
    filename = 'reference_initial_'+fragment+'.'+format
    filename = get_initial_reference_foldername(pname)+filename
    return filename


def get_primers_filename(pname, format='fasta'):
    '''Get the filename with the patient-specific primers'''
    filename = 'primers_'+pname+'.'+format
    filename = get_foldername(pname)+filename
    return filename


def get_mapped_to_initial_foldername(pname, samplename_pat, PCR=1):
    '''Get the folder of mapped reads to initial reference'''
    fn = get_sample_foldername(pname, samplename_pat, PCR=PCR)+'mapped_to_initial/'
    return fn


def get_mapped_to_initial_filename(pname, samplename_pat, 
                                   samplename, fragment, type='bam',
                                   PCR=1,
                                   part=None, unsorted=False,
                                   only_chunk=None):
    '''Get the filename of the mapped reads to initial reference'''
    filename = samplename+'_'+fragment
    if part is not None:
        filename = filename+'_part'+str(part)
    elif unsorted:
        filename = filename+'_unsorted'
    if only_chunk is not None:
        filename = filename+'_chunk_'+str(only_chunk)
    filename = filename+'.'+type
    filename = get_mapped_to_initial_foldername(pname, samplename_pat, PCR=PCR)+filename
    return filename


def get_mapped_filtered_filename(pname, samplename_pat, fragment, type='bam', PCR=1,
                                 decontaminated=False):
    '''Get the filename of the mapped and filtered reads to initial reference'''
    filename = fragment
    if not decontaminated:
        filename = filename+'_to_decontaminate'
    filename = filename+'.'+type
    filename = get_mapped_to_initial_foldername(pname, samplename_pat, PCR=PCR)+filename
    return filename


def get_initial_index_filename(pname, fragment, ext=True):
    '''Get the index filename, with or w/o extension'''
    filename = 'reference_initial_'+fragment
    filename = 'hash/'+filename
    filename = get_initial_reference_foldername(pname)+filename
    if ext:
        filename = filename+'.stidx'
    return filename


def get_initial_hash_filename(pname, fragment, ext=True):
    '''Get the index filename, with or w/o extension'''
    filename = 'reference_initial_'+fragment
    filename = 'hash/'+filename
    filename = get_initial_reference_foldername(pname)+filename
    if ext:
        filename = filename+'.sthash'
    return filename


def get_consensi_alignment_filename(pname, fragment, format='fasta'):
    '''Get the MSA of all consensi of the patient, sorted by time point'''
    filename = 'consensi_alignment_'+fragment+'.'+format
    filename = get_foldername(pname)+filename
    return filename


def get_consensi_tree_filename(pname, fragment):
    '''Get the newick filename of the consensus tree'''
    filename = 'consensi_tree_'+fragment+'.newick'
    filename = get_foldername(pname)+filename
    return filename


def get_consensi_alignment_genomewide_filename(pname):
    '''Get the MSA of all consensi of the patient, sorted by time point'''
    filename = 'consensi_alignment_genomewide.fasta'
    filename = get_foldername(pname)+filename
    return filename


def get_allele_counts_filename(pname, samplename_pat, fragment, PCR=1, qual_min=30):
    '''Get the filename of the allele counts for a patient sample'''
    filename = 'allele_counts_'+fragment+'_qual'+str(qual_min)+'+'+'.npy'
    filename = get_sample_foldername(pname, samplename_pat, PCR=PCR)+filename
    return filename


def get_allele_cocounts_filename(pname, samplename_pat, fragment, PCR=1, qual_min=30,
                                 compressed=True):
    '''Get the matrix of allele cocounts on the initial reference'''
    filename = 'allele_cocounts_'+fragment+'_qual'+str(qual_min)+'+'+'.'
    if compressed:
        filename = filename+'npz'
    else:
        filename = filename+'npy'

    filename = get_sample_foldername(pname, samplename_pat, PCR=PCR)+filename
    return filename


def get_consensus_filename(pname, samplename_pat, fragment, PCR=1, format='fasta'):
    '''Get the filename of the consensus of a patient sample'''
    filename = 'consensus_'+fragment+'.'+format
    filename = get_sample_foldername(pname, samplename_pat, PCR=PCR)+filename
    return filename


def get_allele_count_trajectories_filename(pname, fragment):
    '''Get the matrix with allele counts on the initial reference'''
    filename = 'allele_counts_trajectories_'+fragment+'.npz'
    filename = get_foldername(pname)+filename
    return filename


def get_allele_frequency_trajectories_filename(pname, fragment):
    '''Get the matrix with allele frequencies on the initial reference'''
    filename = 'allele_frequency_trajectories_'+fragment+'.npz'
    filename = get_foldername(pname)+filename
    return filename


def get_coordinate_map_filename(pname, fragment, refname='HXB2'):
    '''Get the filename of the map to HXB2 or other ref coordinates'''
    filename = 'map_coord_to_'+refname+'_'+fragment+'.dat'
    filename = get_foldername(pname)+filename
    return filename


def get_SFS_filename(pnames, fragments, suffix=None):
    '''Get the filename of the SFS'''
    filename = 'SFS'
    if len(fragments) == 6:
        filename = filename+'_allfrags'
    else:
        filename = filename+'_'+'_'.join(fragments)
    
    if len(pnames) == 1:
        pname = pnames[0]
    else:
        pname = 'all'
        filename = filename+'_'+'_'.join(pnames)

    if suffix is not None:
        filename = filename+'_'+suffix

    filename = filename+'.npz'

    filename = get_foldername(pname)+'sfs/'+filename
    return filename


def get_propagator_filename(pnames, fragments, dt, suffix=None):
    '''Get the filename of the propagator'''
    filename = 'propagator'
    if len(fragments) == 6:
        filename = filename+'_allfrags'
    else:
        filename = filename+'_'+'_'.join(fragments)
    
    if len(pnames) == 1:
        pname = pnames[0]
    else:
        pname = 'all'
        filename = filename+'_'+'_'.join(pnames)

    filename = filename+'_dt_'+'_'.join(map(str, dt))

    if suffix is not None:
        filename = filename+'_'+suffix

    filename = filename+'.npz'

    filename = get_foldername(pname)+'propagator/'+filename
    return filename


def get_divergence_trajectories_local_filename(pname, fragment):
    '''Get filename of the trajectories of local divergence'''
    filename = 'divergence_trajectories_local_'+fragment+'.npz'
    filename = get_foldername(pname)+filename
    return filename
    

def get_diversity_trajectories_local_filename(pname, fragment):
    '''Get filename of the trajectories of local diversity'''
    filename = 'diversity_trajectories_local_'+fragment+'.npz'
    filename = get_foldername(pname)+filename
    return filename
    



# FIGURES
def get_figure_folder(pname):
    '''Get the folder for figures for this sample'''
    folder = get_foldername(pname)+'figures/'
    return folder


def get_allele_frequency_trajectory_figure_filename(pname, fragment_or_gene, format='png'):
    '''Get the filename of the plot of allele frequency trajectories'''
    folder = get_figure_folder(pname)
    fn = folder+'allele_freq_traj_'+fragment_or_gene+'.'+format
    return fn


def get_correlation_PCR1_PCR2_aft_figure_filename(pname, fragment, samplename, format='png'):
    '''Get the filename of the figure of allele freq correlation between PCR1 and PCR2'''
    folder = get_figure_folder(pname)
    fn = folder+'correlation_PCR1_PCR2_allele_frequencies_'+fragment+'_'+samplename+'.'+format
    return fn


def get_coverage_to_initial_figure_filename(pname, fragment, format='png'):
    '''Get the filename of the figure of coverage along the infection'''
    folder = get_figure_folder(pname)
    fn = folder+'coverage_'+fragment+'.'+format
    return fn


def get_tree_consensi_figure_filename(pname, fragment, format='png'):
    '''Get the filename of the figure of the tree of consensi'''
    folder = get_figure_folder(pname)
    fn = folder+'tree_consensi_'+fragment+'.'+format
    return fn


def get_crosscontamination_figure_filename(fragment, format='pdf'):
    '''Get the filename of the figure with the cross-contamination matrix'''
    folder = get_figure_folder('all')
    fn = folder+'crosscontamination_'+fragment+'.'+format
    return fn


def get_physiological_figure_filename(pname, format='png'):
    '''Get the filename of the figure with the viral load and CD4+ counts'''
    folder = get_figure_folder(pname)
    fn = folder+'physiological.'+format
    return fn


def get_SFS_figure_filename(pname, fragments, format='png'):
    '''Get the filename of the figure of the site frequency spectrum'''
    folder = get_figure_folder(pname)
    fn = folder+'SFS_'+'-'.join(fragments)+'.'+format
    return fn



# SUMMARY
def get_map_initial_summary_filename(pname, samplename_pat, samplename, fragment, PCR=1, only_chunk=None):
    '''Get the filename of the summary of the map to initial reference'''
    filename = 'summary_map_initial_'+samplename+'_'+fragment
    if only_chunk is not None:
        filename = filename+'_chunk_'+str(only_chunk)
    filename = filename+'.txt'
    filename = get_sample_foldername(pname, samplename_pat, PCR=PCR)+'mapped_to_initial/'+filename
    return filename


def get_paste_mapped_chunks_initial_summary_filename(pname, samplename_pat, samplename, fragment, PCR=1):
    '''Get the filename of the summary of the pasting of chunk after mapping to initial reference'''
    filename = 'summary_paste_mapped_initial_'+samplename+'_'+fragment+'.txt'
    filename = get_sample_foldername(pname, samplename_pat, PCR=PCR)+'mapped_to_initial/'+filename
    return filename


def get_filter_mapped_init_summary_filename(pname, samplename_pat, fragment, PCR=1):
    '''Get the filename of the summary of the post-map filtering'''
    filename= 'summary_filter_mapped_init_'+fragment+'.txt'
    filename = get_sample_foldername(pname, samplename_pat, PCR=PCR)+'mapped_to_initial/'+filename
    return filename


def get_decontaminate_summary_filename(pname, samplename_pat, fragment, PCR=1):
    '''Get the filename of the summary of the decontamination'''
    filename = 'summary_decontaminate_'+fragment+'.txt'
    filename = get_sample_foldername(pname, samplename_pat, PCR=PCR)+'mapped_to_initial/'+filename
    return filename
