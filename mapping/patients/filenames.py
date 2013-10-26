# vim: fdm=marker
'''
author:     Fabio Zanini
date:       25/10/13
content:    Module containing all filenames for the patient analysis in one place.
'''
# Modules
import os



# Functions
def get_foldername(pname, root_data_folder=root_data_folder):
    '''Get the folder name of the data from a patient'''
    foldername = 'patients/'+patient+'/'
    foldername = root_data_folder+foldername
    return foldername


def get_initial_consensus_filename(pname, fragment,
                                   root_data_folder=root_data_folder):
    '''Get the filename of the initial consensus for a patient'''
    filename = 'consensus_initial_'+fragment+'.fasta'
    filename = get_foldername(pname, root_data_folder=root_data_folder)+filename
    return filename


def get_mapped_to_initial_filename(pname, sample, fragment, type='bam',
                                   part=None, unsorted=False,
                                   root_data_folder=root_data_folder):
    '''Get the filename of the mapped reads to initial consensus'''
    filename = fragment
    if part is not None:
        filename = filename+'_part'+str(part)
    elif unsorted:
        filename = filename+'_unsorted'
    filename = filename+'.'+type
    filename = sample+'/mapped_to_initial/'+filename
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


def get_consensi_alignment_genomewide_filename(pname,
                                               root_data_folder=root_data_folder):
    '''Get the MSA of all consensi of the patient, sorted by time point'''
    filename = 'consensi_alignment_genomewide.fasta'
    filename = get_foldername(pname, root_data_folder=root_data_folder)+filename
    return filename
