#!/ebio/ag-neher/share/programs/EPD/bin/python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/08/13
content:    Map the reads to their own consensus, produced with an iterative
            mapping of HIV onto itself (HXB2 at first). This is the final mapping.
'''
# Modules
import os
import sys
import argparse
import numpy as np
import subprocess as sp
from adapter_info import load_adapter_table, foldername_adapter
from mapping.mapping_utils import stampy_bin, subsrate, bwa_bin, convert_sam_to_bam
from mapping.filenames import get_consensus_filename, get_mapped_filename,\
        get_read_filenames



# Globals
VERBOSE = 1
# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']


# Cluster submit
import mapping
JOBDIR = mapping.__path__[0].rstrip('/')+'/'
JOBSCRIPT = JOBDIR+'map_to_consensus.py'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
# Different times based on subsample flag
cluster_time = ['23:59:59', '0:59:59']
vmem = '8G'


# Functions
def get_index_file(data_folder, adaID, fragment, subsample=False, ext=True):
    '''Get the index filename, with or w/o extension'''
    filename = 'consensus_'+fragment
    filename = 'hash/'+filename
    filename = foldername_adapter(adaID)+filename
    if subsample:
        filename = 'subsample/'+filename
    if ext:
        filename = filename+'.stidx'
    return data_folder+filename


def get_hash_file(data_folder, adaID, fragment, subsample=False, ext=True):
    '''Get the index filename, with or w/o extension'''
    filename = 'consensus_'+fragment
    filename = 'hash/'+filename
    filename = foldername_adapter(adaID)+filename
    if subsample:
        filename = 'subsample/'+filename
    if ext:
        filename = filename+'.sthash'
    return data_folder+filename


def make_index_and_hash(data_folder, adaID, fragment, subsample=False, VERBOSE=0):
    '''Make index and hash files for consensus'''
    # Make folder if necessary
    dirname =  os.path.dirname(get_hash_file(data_folder, adaID, 'F0',
                                             subsample=subsample))
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    # 1. Make genome index file
    if not os.path.isfile(get_index_file(data_folder, adaID, fragment,
                                         subsample=subsample, ext=True)):
        if VERBOSE:
            print 'Build index: '+'{:02d}'.format(adaID)+' '+fragment
        sp.call([stampy_bin,
                 '--species="HIV fragment '+fragment+'"',
                 '-G', get_index_file(data_folder, adaID, fragment,
                                      subsample=subsample, ext=False),
                 get_consensus_filename(data_folder, adaID, fragment,
                                        subsample=subsample),
                 ])
    
    # 2. Build a hash file
    if not os.path.isfile(get_hash_file(data_folder, adaID, fragment,
                                        subsample=subsample, ext=True)):
        if VERBOSE:
            print 'Build hash: '+'{:02d}'.format(adaID)+' '+fragment
        sp.call([stampy_bin,
                 '-g', get_index_file(data_folder, adaID, fragment,
                                      subsample=subsample, ext=False),
                 '-H', get_hash_file(data_folder, adaID, fragment,
                                     subsample=subsample, ext=False),
                 ])


def make_bwa_hash(data_folder, adaID, fragment, subsample=False, VERBOSE=0):
    '''Make hash tables for BWA'''
    # Make folder if necessary
    dirname =  os.path.dirname(get_hash_file(data_folder, adaID, 'F0',
                                             subsample=subsample))
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    # Call bwa index
    cons_filename = get_consensus_filename(data_folder, adaID, fragment,
                                           subsample=subsample)
    if VERBOSE:
        print 'Build BWA index: '+'{:02d}'.format(adaID)+' '+fragment
    sp.call([bwa_bin,
             'index',
             cons_filename,
             ])

    # Move the hashes into subfolder
    from glob import glob
    import shutil
    shutil.copy(cons_filename, dirname)
    hash_files = glob(cons_filename+'.*')
    for hash_file in hash_files:
        shutil.copy(hash_file, dirname)
        os.remove(hash_file)


def map_bwa(data_folder, adaID, fragment, subsample=False, VERBOSE=0):
    '''Map using BWA'''
    index_prefix = get_hash_file(data_folder, adaID, fragment,
                                 subsample=subsample, ext=False)+'.fasta'

    if VERBOSE:
        print 'Map via BWA: '+'{:02d}'.format(adaID)+' '+fragment
    mapped_filename = get_mapped_filename(data_folder, adaID, fragment,
                                          type='sam', subsample=subsample,
                                          bwa=True)
    # Make folder if necessary
    dirname = os.path.dirname(mapped_filename)
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    # Map
    with open(mapped_filename, 'w') as f:
        call_list = [bwa_bin,
                     'mem',
                     index_prefix] +\
                    get_read_filenames(data_folder, adaID, subsample=subsample,
                                       premapped=True, fragment=fragment)
        if VERBOSE >=2:
            print ' '.join(call_list)
        sp.call(call_list, stdout=f)


def map_stampy_after_bwa(data_folder, adaID, fragment, subsample=False, VERBOSE=0):
    '''Map using stampy after BWA'''
    if VERBOSE:
        print 'Map via stampy after BWA: '+'{:02d}'.format(adaID)+' '+fragment

    bwa_filename = get_mapped_filename(data_folder, adaID, fragment,
                                       type='bam', subsample=subsample,
                                       bwa=True)

    if not os.path.isfile(bwa_filename):
        if VERBOSE >= 2:
            print 'Converting SAM to BAM'
        convert_sam_to_bam(bwa_filename)

    mapped_filename = get_mapped_filename(data_folder, adaID, fragment,
                                          type='sam', subsample=subsample)
    # Map
    call_list = [stampy_bin,
                 '-g', get_index_file(data_folder, adaID, fragment,
                                      subsample=subsample, ext=False),
                 '-h', get_hash_file(data_folder, adaID, fragment,
                                     subsample=subsample, ext=False), 
                 '-o', mapped_filename,
                 '--substitutionrate='+subsrate,
                 '--bamkeepgoodreads',
                 '-M', bwa_filename]
    if VERBOSE >=2:
        print ' '.join(call_list)
    sp.call(call_list)

             

def fork_self(data_folder, adaID, fragment, subsample=False, VERBOSE=3):
    '''Fork self for each adapter ID and fragment'''
    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'map '+'{:02d}'.format(adaID)+' '+fragment,
                 '-l', 'h_rt='+cluster_time[subsample],
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--adaIDs', adaID,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                ]
    if subsample:
        qsub_list.append('--subsample')
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    sp.call(qsub_list)

    

# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Divide HIV reads into fragments')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--subsample', action='store_true',
                        help='Apply only to a subsample of the reads')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')

    args = parser.parse_args()
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    subsample = args.subsample
    submit = args.submit

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = load_adapter_table(data_folder)['ID']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Iterate over all requested samples
    for adaID in adaIDs:
        for fragment in fragments:

            # Submit to the cluster self if requested
            if submit:
                fork_self(data_folder, adaID, fragment,
                          subsample=subsample, VERBOSE=VERBOSE)
                continue

            # Make BWA hashes
            make_bwa_hash(data_folder, adaID, fragment,
                          subsample=subsample, VERBOSE=VERBOSE)
    
            # Map via BWA first
            map_bwa(data_folder, adaID, fragment,
                    subsample=subsample, VERBOSE=VERBOSE)
    
            # Make stampy hashes
            make_index_and_hash(data_folder, adaID, fragment,
                                subsample=subsample, VERBOSE=VERBOSE)
            
            # Map via stampy afterwards
            map_stampy_after_bwa(data_folder, adaID, fragment,
                                 subsample=subsample, VERBOSE=VERBOSE)
    
