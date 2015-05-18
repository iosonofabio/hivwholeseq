# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/11/14
content:    Module containing filenames for the website.
'''
# Modules
from hivwholeseq.filenames import root_data_folder


# Globals
root_website_folder = root_data_folder+'website/'



# Functions
def get_foldername(name):
    '''Get the folder name of an observable'''
    foldername = root_website_folder+name+'/'
    return foldername


def get_sample_table_filename(pname, kind):
    '''Get filename of the same table'''
    filename = 'samples_'+kind+'_'+pname+'.tsv'
    filename = get_foldername('tables')+filename
    return filename


def get_tree_figure_filename(pname, region, format='svg'):
    '''Get filename of the example tree for the patient page'''
    filename = 'tree_'+pname+'_'+region+'_example.'+format
    filename = get_foldername('trees')+filename
    return filename


def get_consensi_tree_filename(pname, fragment, format='json'):
    '''Get the filename of a consensus tree'''
    filename = 'consensi_tree_'+pname+'_'+fragment+'.'+format
    filename = get_foldername('trees')+filename
    return filename


def get_consensi_alignment_filename(pname, fragment, format='fasta'):
    '''Get the filename of a consensus alignment'''
    filename = 'consensi_alignment_'+pname+'_'+fragment+'.'+format
    filename = get_foldername('alignments')+filename
    return filename


def get_haplotype_tree_filename(pname, fragment, format='json'):
    '''Get the filename of a haplotype tree'''
    filename = 'haplotype_tree_'+pname+'_'+fragment+'.'+format
    filename = get_foldername('trees')+filename
    return filename


def get_haplotype_alignment_filename(pname, fragment, format='fasta'):
    '''Get the filename of a haplotype alignment'''
    filename = 'haplotype_alignment_'+pname+'_'+fragment+'.'+format
    filename = get_foldername('alignments')+filename
    return filename


def get_viral_load_filename(pname):
    '''Get the filename of the viral load'''
    filename = 'viral_load_'+pname+'.dat'
    filename = get_foldername('physiological')+filename
    return filename


def get_cell_count_filename(pname):
    '''Get the filename of the CD4+ counts'''
    filename = 'cell_count_'+pname+'.dat'
    filename = get_foldername('physiological')+filename
    return filename


def get_divergence_filename(pname, fragment):
    '''Get the filename of genetic divergence'''
    filename = 'divergence_'+pname+'_'+fragment+'.dat'
    filename = get_foldername('divdiv')+filename
    return filename


def get_diversity_filename(pname, fragment):
    '''Get the filename of genetic diversity'''
    filename = 'diversity_'+pname+'_'+fragment+'.dat'
    filename = get_foldername('divdiv')+filename
    return filename


def get_coverage_filename(pname, fragment):
    '''Get the filename of coverage'''
    filename = 'coverage_'+pname+'_'+fragment+'.npz'
    filename = get_foldername('coverage')+filename
    return filename


def get_allele_count_trajectories_filename(pname, fragment, format='npz'):
    '''Get the filename of allele count trajectories'''
    filename = 'allele_counts_'+pname+'_'+fragment+'.'+format
    filename = get_foldername('single_nucleotide_variants')+filename
    return filename


def get_patient_reference_filename(pname, fragment='genomewide', format='gb'):
    '''Get the filename of the patient reference'''
    filename = 'reference_'+pname+'_'+fragment+'.'+format
    filename = get_foldername('sequences')+filename
    return filename
    

def get_coordinate_map_filename(pname, fragment='genomewide', refname='HXB2',
                                format='tsv'):
    '''Get the filename of the patient reference'''
    filename = 'coordinate_map_'+pname+'_'+refname+'_'+fragment+'.'+format
    filename = get_foldername('coordinate_maps')+filename
    return filename


def get_SFS_filename(pnames, fragments, suffix=None):
    '''Get the filename of the SFS'''
    filename = 'sfs'
    if len(fragments) == 6:
        filename = filename+'_allfrags'
    else:
        filename = filename+'_'+'_'.join(fragments)
    
    filename = filename+'_'+'_'.join(pnames)

    if suffix is not None:
        filename = filename+'_'+suffix

    filename = filename+'.npz'

    filename = get_foldername('one_site')+filename
    return filename


def get_propagator_filename(pnames, fragments, dt, suffix=None):
    '''Get the filename of the propagator'''
    filename = 'propagator'
    if len(fragments) == 6:
        filename = filename+'_allfrags'
    else:
        filename = filename+'_'+'_'.join(fragments)
    
    filename = filename+'_'+'_'.join(pnames)

    filename = filename+'_dt_'+'_'.join(map(str, dt))

    if suffix is not None:
        filename = filename+'_'+suffix

    filename = filename+'.npz'

    filename = get_foldername('one_site')+filename
    return filename


def get_divergence_trajectories_filename(pname, fragment):
    '''Get the filename of the local divergence trajectories'''
    filename = 'divergence_trajectory_'+pname+'_'+fragment+'.npz'
    filename = get_foldername('divdiv')+filename
    return filename


def get_diversity_trajectories_filename(pname, fragment):
    '''Get the filename of the local diversity trajectories'''
    filename = 'diversity_trajectory_'+pname+'_'+fragment+'.npz'
    filename = get_foldername('divdiv')+filename
    return filename


def get_reads_filename(pname, fragment, it, format='bam'):
    '''Get the filename of the mapped reads'''
    filename = fragment+'.'+format
    filename = str(it+1)+'/'+filename
    filename = pname+'/samples/'+filename
    filename = get_foldername('patients')+filename
    return filename


def get_timeline_filename(pname):
    '''Get the filename of the patient timeline'''
    filename = 'timeline.tsv'
    filename = get_foldername('patients')+pname+'/'+filename
    return filename


def get_ntemplates_filename(pname):
    '''Get the filename of the template numbers'''
    filename = 'ntemplates_'+pname+'.dat'
    filename = get_foldername('physiological')+filename
    return filename
