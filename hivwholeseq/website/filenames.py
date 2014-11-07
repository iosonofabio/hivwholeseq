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


def get_consensi_tree_filename(pname, fragment, format='newick'):
    '''Get the newick filename of a consensus tree'''
    filename = 'consensi_tree_'+pname+'_'+fragment+'.'+format
    filename = get_foldername('trees')+filename
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
    '''Get the filename of genetic diversity'''
    filename = 'coverage_'+pname+'_'+fragment+'.npz'
    filename = get_foldername('one_site')+filename
    return filename

