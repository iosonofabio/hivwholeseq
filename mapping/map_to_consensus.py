#!/usr/bin/env python
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
import time
import argparse
import pysam
import numpy as np
import subprocess as sp

from mapping.datasets import MiSeq_runs
from mapping.adapter_info import load_adapter_table, foldername_adapter
from mapping.mapping_utils import stampy_bin, subsrate, bwa_bin, convert_sam_to_bam, \
        convert_bam_to_sam
from mapping.filenames import get_consensus_filename, get_mapped_filename,\
        get_read_filenames



# Globals
# Stampy parameters
stampy_gapopen = 60	        # Default: 40
stampy_gapextend = 5 	    # Default: 3
stampy_sensitive = True    # Default: False


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
def fork_self(miseq_run, adaID, fragment, subsample=False, VERBOSE=3, bwa=False,
              threads=1):
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
                 '--run', miseq_run,
                 '--adaIDs', adaID,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                 '--threads', threads,
                ]
    if subsample:
        qsub_list.append('--subsample')
    if bwa:
        qsub_list.append('--bwa')
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


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


def make_output_folders(data_folder, adaID, subsample=False, VERBOSE=0, bwa=False):
    '''Make the output folders if necessary for hash and map'''
    hash_foldername = os.path.dirname(get_hash_file(data_folder, adaID, 'F0',
                                                    subsample=subsample))
    map_foldername = os.path.dirname(get_mapped_filename(data_folder, adaID, 'F0',
                                     subsample=subsample))
    foldernames = [hash_foldername, map_foldername]

    # Add BWA folder if requested
    if bwa:
        bwa_foldername = os.path.dirname(get_mapped_filename(data_folder, adaID, 'F0',
                                         subsample=subsample, bwa=True))
        foldernames.append(bwa_foldername)

    # Make the folders
    for dirname in foldernames:
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        if VERBOSE:
            print 'Folder created:', dirname


def make_index_and_hash(data_folder, adaID, fragment, subsample=False, VERBOSE=0):
    '''Make index and hash files for consensus'''
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
    output_filename = get_mapped_filename(data_folder, adaID, fragment,
                                          type='sam', subsample=subsample,
                                          bwa=True)

    # Map
    with open(output_filename, 'w') as f:
        call_list = [bwa_bin,
                     'mem',
                     index_prefix] +\
                    get_read_filenames(data_folder, adaID, subsample=subsample,
                                       premapped=True, fragment=fragment)
        if VERBOSE >=2:
            print ' '.join(call_list)
        sp.call(call_list, stdout=f)


def map_stampy(data_folder, adaID, fragment, subsample=False, VERBOSE=0, bwa=False,
               threads=1):
    '''Map using stampy, either directly or after BWA'''
    # Get input filenames
    if bwa:
        if VERBOSE:
            print 'Map via stampy after BWA: '+'{:02d}'.format(adaID)+' '+fragment

        bwa_filename = get_mapped_filename(data_folder, adaID, fragment,
                                           type='bam', subsample=subsample,
                                           bwa=True)
        if not os.path.isfile(bwa_filename):
            if VERBOSE >= 2:
                print 'Converting SAM to BAM'
            convert_sam_to_bam(bwa_filename)
        readfiles = [bwa_filename]

    else:
        if VERBOSE:
            print 'Map via stampy (no BWA): '+'{:02d}'.format(adaID)+' '+fragment
        readfiles = get_read_filenames(data_folder, adaID, fragment=fragment,
                                       subsample=subsample, premapped=True)

    # parallelize if requested
    if threads == 1:

        # Get output filename
        output_filename = get_mapped_filename(data_folder, adaID, fragment,
                                              type='sam', subsample=subsample)

        # Map
        call_list = [stampy_bin,
                     '-g', get_index_file(data_folder, adaID, fragment,
                                          subsample=subsample, ext=False),
                     '-h', get_hash_file(data_folder, adaID, fragment,
                                         subsample=subsample, ext=False), 
                     '-o', output_filename,
                     '--substitutionrate='+subsrate,
                     '--gapopen', stampy_gapopen,
                     '--gapextend', stampy_gapextend]
        if bwa:
            call_list.append('--bamkeepgoodreads')
        if stampy_sensitive:
            call_list.append('--sensitive')
        call_list = call_list + ['-M'] + readfiles
        call_list = map(str, call_list)
        if VERBOSE >=2:
            print ' '.join(call_list)
        sp.call(call_list)

    else:

        # Submit map script
        jobs_done = np.zeros(threads, bool)
        job_IDs = np.zeros(threads, 'S30')
        for j in xrange(threads):
    
            # Get output filename
            output_filename =  get_mapped_filename(data_folder, adaID, fragment,
                                               type='sam', subsample=subsample,
                                               part=(j+1))
            # Map
            call_list = ['qsub','-cwd',
                         '-b', 'y',
                         '-S', '/bin/bash',
                         '-o', JOBLOGOUT,
                         '-e', JOBLOGERR,
                         '-N', 'm '+'{:02d}'.format(adaID)+fragment+' p'+str(j+1),
                         '-l', 'h_rt='+cluster_time[subsample],
                         '-l', 'h_vmem='+vmem,
                         stampy_bin,
                         '-g', get_index_file(data_folder, adaID, fragment,
                                              subsample=subsample, ext=False),
                         '-h', get_hash_file(data_folder, adaID, fragment,
                                             subsample=subsample, ext=False), 
                         '-o', output_filename,
                         '--processpart='+str(j+1)+'/'+str(threads),
                         '--substitutionrate='+subsrate,
                         '--gapopen', stampy_gapopen,
                         '--gapextend', stampy_gapextend]
            if bwa:
                call_list.append('--bamkeepgoodreads')
            if stampy_sensitive:
                call_list.append('--sensitive')
            call_list = call_list + ['-M'] + readfiles
            call_list = map(str, call_list)
            if VERBOSE >= 2:
                print ' '.join(call_list)
            job_ID = sp.check_output(call_list)
            job_ID = job_ID.split()[2]
            job_IDs[j] = job_ID

        # Check output
        while not jobs_done.all():

            # Sleep 10 secs
            time.sleep(10)

            # Get the output of qstat to check the status of jobs
            qstat_output = sp.check_output(['qstat'])
            qstat_output = qstat_output.split('\n')[:-1] # The last is an empty line
            if len(qstat_output) < 3:
                jobs_done[:] = True
                break
            else:
                qstat_output = [line.split()[0] for line in qstat_output[2:]]

            for j in xrange(threads):
                if jobs_done[j]:
                    continue

                if job_IDs[j] not in qstat_output:
                    jobs_done[j] = True

        # Merge output files
        output_file_parts = [get_mapped_filename(data_folder, adaID, fragment,
                                             type='bam', subsample=subsample,
                                             part=(j+1))
                              for j in xrange(threads)]
        for output_file in output_file_parts:
            convert_sam_to_bam(output_file)
        output_filename = get_mapped_filename(data_folder, adaID, fragment, type='bam',
                                              subsample=subsample, unsorted=True)
        pysam.merge(output_filename, *output_file_parts)

        # Sort the file by read names (to ensure the pair_generator)
        output_filename_sorted = get_mapped_filename(data_folder, adaID, fragment,
                                                     type='bam',
                                                     subsample=subsample,
                                                     unsorted=False)

        # Note: we exclude the extension and the option -f because of a bug in samtools
        pysam.sort('-n', output_filename, output_filename_sorted[:-4])

        # Make SAM out of the BAM for checking
        output_filename = get_mapped_filename(data_folder, adaID, fragment,
                                          type='sam', subsample=subsample)
        convert_bam_to_sam(output_filename)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Divide HIV reads into fragments')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
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
    parser.add_argument('--bwa', action='store_true',
                        help='Use BWA for premapping?')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads to use for mapping')

    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    subsample = args.subsample
    submit = args.submit
    bwa = args.bwa
    threads = args.threads

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

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

        # Make output folders if necessary
        make_output_folders(data_folder, adaID, subsample=subsample)

        for fragment in fragments:

            # Submit to the cluster self if requested
            if submit:
                fork_self(miseq_run, adaID, fragment,
                          subsample=subsample, VERBOSE=VERBOSE, bwa=bwa,
                          threads=threads)
                continue

            if bwa:
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
            map_stampy(data_folder, adaID, fragment,
                       subsample=subsample, VERBOSE=VERBOSE, bwa=bwa,
                       threads=threads)
    
