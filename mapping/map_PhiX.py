#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       03/08/13
content:    Map the unclassified reads to phiX as a control.
'''
# Modules
import os
import subprocess as sp
import argparse

from mapping.datasets import MiSeq_runs
from mapping.filenames import get_mapped_phix_filename, get_phix_filename, \
        get_unclassified_reads_filenames
from mapping.mapping_utils import stampy_bin


# Globals
# Submit vars
JOBSUBMIT = os.path.realpath(__file__)      # This file
JOBDIR = os.path.dirname(JOBSUBMIT)+'/'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'

# Cluster parameters
cluster_time = '23:59:59'
vmem = '8G'
    

# Functions
def get_phix_index_filename(data_folder, ext=True):
    '''Get index file for phiX'''
    dirname = os.path.dirname(get_mapped_phix_filename(data_folder)).rstrip('/')+'/'
    filename = dirname+'hash/phix'
    if ext:
        filename = filename+'.stidx'
    return filename


def get_phix_hash_filename(data_folder, ext=True):
    '''Get hash file for phiX'''
    dirname = os.path.dirname(get_mapped_phix_filename(data_folder)).rstrip('/')+'/'
    filename = dirname+'hash/phix'
    if ext:
        filename = filename+'.sthash'
    return filename


def make_output_folders(data_folder, VERBOSE=0):
    '''Make the phix folder if necessary'''
    for dirname in map(os.path.dirname,
                       [get_mapped_phix_filename(data_folder),
                        get_phix_index_filename(data_folder),
                        get_phix_hash_filename(data_folder)]):

        if not os.path.isdir(dirname):
            os.mkdir(dirname)
            if VERBOSE >= 2:
                print 'Folder created:', dirname



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Map reads to PhiX')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')

    args = parser.parse_args()
    miseq_run = args.run
    VERBOSE = args.verbose
    submit = args.submit

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # Make output folder
    make_output_folders(data_folder)

    # Get reference sequence filename
    ref_filename = get_phix_filename()

    # 1. Make genome index file
    if not os.path.isfile(get_phix_index_filename(data_folder)):
        sp.call([stampy_bin,
                 '--species=phiX',
                 '-G', get_phix_index_filename(data_folder, ext=False),
                 ref_filename,
                 ])
    
    # 2. Build a hash file
    if not os.path.isfile(get_phix_hash_filename(data_folder)):
        sp.call([stampy_bin,
                 '-g', get_phix_index_filename(data_folder, ext=False),
                 '-H', get_phix_hash_filename(data_folder, ext=False),
                 ])
    
    # 3. Map using stampy
    input_filenames = get_unclassified_reads_filenames(data_folder, filtered=True)[:2]
    call_list = [stampy_bin,
                 '-g', get_phix_index_filename(data_folder, ext=False),
                 '-h', get_phix_hash_filename(data_folder, ext=False), 
                 '-o', get_mapped_phix_filename(data_folder, type='sam'),
                 '-M'] + input_filenames

    if submit:
        call_list = ['qsub','-cwd',
                     '-o',JOBLOGOUT,
                     '-e',JOBLOGERR,
                     '-N', 'map phiX',
                     '-l', 'h_rt='+cluster_time,
                     '-l', 'h_vmem='+vmem] + call_list

    call_list = map(str, call_list)
    if VERBOSE:
        print ' '.join(call_list)
    sp.call(call_list)
