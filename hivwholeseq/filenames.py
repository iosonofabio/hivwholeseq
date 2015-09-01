# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/10/14
content:    Module containing basic paths for the whole analysis.
'''
# Modules
import os
from warnings import warn

import hivwholeseq as self
from .utils.generic import which


# Globals
root_data_folder = 'HIVWHOLESEQ_ROOT_DATA_FOLDER'
if root_data_folder not in os.environ:
    raise ValueError('Environment variable not found: '+root_data_folder)
elif not os.path.isdir(os.environ[root_data_folder]):
    raise IOError('Root data folder is not a folder')
else:
    root_data_folder = os.environ[root_data_folder].rstrip('/')+'/'

stampy_bin = 'STAMPY_BIN'
if stampy_bin in os.environ:
    if not os.path.isfile(os.environ[stampy_bin]):
        raise IOError('Stampy bin is not a file')
    else:
        stampy_bin = os.environ[stampy_bin]
else:
    _locs = which('stampy')
    if len(_locs):
        stampy_bin = _locs[0]
    else:
        raise ValueError('Stampy not found. Install it in your PATH or set '+
                         'the environment variable '+stampy_bin+'.')


fasttree_bin = 'FASTTREE_BIN'
if fasttree_bin in os.environ:
    if not os.path.isfile(os.environ[fasttree_bin]):
        raise IOError('FastTree bin is not a file')
    else:
        fasttree_bin = os.environ[fasttree_bin]
else:
    _locs = which('FastTree')
    if len(_locs):
        fasttree_bin = _locs[0]
    else:
        warn('FastTree not found. Install it in your PATH or set '+
             'the environment variable '+fasttree_bin+'.')
 

tmp_folder = root_data_folder+'tmp/'
reference_folder = root_data_folder+'reference/'
theory_folder = root_data_folder+'theory/'
table_folder = self.__path__[0] + '/data/'
table_filename = table_folder+'HIV_reservoir_all.xlsx'



# Functions
def get_custom_reference_filename(reference, format='fasta'):
    '''Get the filename of a custom reference sequence'''
    filename = reference
    filename = filename+'.'+format
    return reference_folder+filename


def get_custom_alignment_filename(aliname, format='fasta'):
    '''Get the filename of a custom alignment'''
    filename = 'alignments/'+aliname
    filename = filename+'.'+format
    return reference_folder+filename
