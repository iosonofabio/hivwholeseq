#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       03/08/13
content:    Map the unclassified reads to phiX as a control.
'''
# Modules
import os
import sys
import argparse
import numpy as np
import pysam

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.sequencing.filenames import get_mapped_phix_filename, get_phix_filename, \
        get_unclassified_reads_filenames
from hivwholeseq.utils.mapping import stampy_bin
from hivwholeseq.utils.mapping import convert_sam_to_bam, convert_bam_to_sam


# Globals
subsrate = '1e-3'

# Cluster submit
import hivwholeseq
JOBDIR = hivwholeseq.__path__[0].rstrip('/')+'/'
JOBLOGERR = JOBDIR+'logerr/'
JOBLOGOUT = JOBDIR+'logout/'
JOBSCRIPT = JOBDIR+'phix/map_phix.py'
cluster_time = '23:59:59'
vmem = '8G'

    

# Functions
def fork_self(seq_run, threads=1, VERBOSE=0):
    '''Fork self'''
    import subprocess as sp

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'mpX'+seq_run,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', seq_run,
                 '--verbose', VERBOSE,
                 '--threads', threads,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


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


def make_index_and_hash(data_folder, VERBOSE=0):
    '''Make index and hash for PhiX'''
    import subprocess as sp

    # 1. Make genome index file
    if not os.path.isfile(get_phix_index_filename(data_folder)):
        sp.call([stampy_bin,
                 '--overwrite',
                 '--species=phiX',
                 '-G', get_phix_index_filename(data_folder, ext=False),
                 get_phix_filename(),
                 ])
    
    # 2. Build a hash file
    if not os.path.isfile(get_phix_hash_filename(data_folder)):
        sp.call([stampy_bin,
                 '--overwrite',
                 '-g', get_phix_index_filename(data_folder, ext=False),
                 '-H', get_phix_hash_filename(data_folder, ext=False),
                 ])


def map_stampy(data_folder, threads=1, VERBOSE=0, subsample=False):
    '''Map using stampy'''
    import subprocess as sp

    input_filenames = get_unclassified_reads_filenames(data_folder)[:2]
    if subsample:
        input_filenames = [r.replace('.fastq', '_subsample.fastq') for r in input_filenames]

    # parallelize if requested
    if threads == 1:
    
        # Map
        call_list = [stampy_bin,
                     '--overwrite',
                     '--sensitive',
                     '-g', get_phix_index_filename(data_folder, ext=False),
                     '-h', get_phix_hash_filename(data_folder, ext=False), 
                     '-o', get_mapped_phix_filename(data_folder, type='sam'),
                     '-M'] + input_filenames

        call_list = map(str, call_list)
        if VERBOSE:
            print ' '.join(call_list)
        sp.call(call_list)

        convert_sam_to_bam(get_mapped_phix_filename(data_folder, type='bam'))
        return

    # Multithreading works as follows: call qsub + stampy, monitor the process
    # IDs with qstat at regular intervals, and finally merge results with pysam
    # Submit map script
    jobs_done = np.zeros(threads, bool)
    job_IDs = np.zeros(threads, 'S30')
    
    # Submit map call
    import time
    import hivwholeseq
    JOBDIR = hivwholeseq.__path__[0].rstrip('/')+'/'
    JOBLOGOUT = JOBDIR+'logout'
    JOBLOGERR = JOBDIR+'logerr'
    cluster_time = ['23:59:59', '1:59:59']
    vmem = '8G'
    for j in xrange(threads):
        call_list = ['qsub','-cwd',
                     '-b', 'y',
                     '-S', '/bin/bash',
                     '-o', JOBLOGOUT,
                     '-e', JOBLOGERR,
                     '-N', 'mpx p'+str(j+1),
                     '-l', 'h_rt='+cluster_time[threads >= 30],
                     '-l', 'h_vmem='+vmem,
                     stampy_bin,
                     '--overwrite',
                     '-g', get_phix_index_filename(data_folder, ext=False),
                     '-h', get_phix_hash_filename(data_folder, ext=False), 
                     '-o', get_mapped_phix_filename(data_folder, type='sam', part=(j+1)),
                     '--processpart='+str(j+1)+'/'+str(threads),
                     '--substitutionrate='+subsrate,
                     '-M'] + input_filenames
        call_list = map(str, call_list)
        if VERBOSE >= 2:
            print ' '.join(call_list)
        job_ID = sp.check_output(call_list)
        job_ID = job_ID.split()[2]
        job_IDs[j] = job_ID

    # Monitor output
    output_file_parts = [get_mapped_phix_filename(data_folder, type='bam', part=(j+1))
                         for j in xrange(threads)]
    time_wait = 5 # secs
    while not jobs_done.all():

        # Sleep some time
        time.sleep(time_wait)

        # Get the output of qstat to check the status of jobs
        qstat_output = sp.check_output(['qstat'])
        qstat_output = qstat_output.split('\n')[:-1] # The last is an empty line
        if len(qstat_output) < 3:
            jobs_done[:] = True
            break
        else:
            qstat_output = [line.split()[0] for line in qstat_output[2:]]

        time_wait = 10 # secs
        for j in xrange(threads):
            if jobs_done[j]:
                continue

            if job_IDs[j] not in qstat_output:
                # Convert to BAM for merging
                if VERBOSE >= 1:
                    print 'Convert premapped reads to BAM for merging, part '+\
                           str(j+1)+' of '+str(threads)
                convert_sam_to_bam(output_file_parts[j])
                # We do not need to wait if we did the conversion (it takes
                # longer than some secs)
                time_wait = 0
                jobs_done[j] = True

    # Concatenate output files
    if VERBOSE >= 1:
        print 'Concatenate premapped reads',
    output_filename = get_mapped_phix_filename(data_folder, type='bam', unsorted=True)
    pysam.cat('-o', output_filename, *output_file_parts)
    if VERBOSE >= 1:
        print 'done.'

    # Sort the file by read names (to ensure the pair_generator)
    # NOTE: we exclude the extension and the option -f because of a bug in samtools
    if VERBOSE >= 1:
        print 'Sort premapped reads'
    output_filename_sorted = get_mapped_phix_filename(data_folder, type='bam', unsorted=False)
    pysam.sort('-n', output_filename, output_filename_sorted[:-4])

    # Reheader the file without BAM -> SAM -> BAM
    if VERBOSE >= 1:
        print 'Reheader premapped reads'
    header_filename = get_mapped_phix_filename(data_folder, type='sam', part=1)
    pysam.reheader(header_filename, output_filename_sorted)



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Map reads to PhiX')
    parser.add_argument('--run', required=True,
                        help='MiSeq run to analyze (e.g. Tue28)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads to use for mapping')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--subsample', action='store_true',
                        help='Execute the script on a subsample of reads only')


    args = parser.parse_args()
    seq_run = args.run
    VERBOSE = args.verbose
    threads = args.threads
    submit = args.submit
    subsample = args.subsample

    # If submit, outsource to the cluster
    if submit:
        fork_self(seq_run, threads=threads, VERBOSE=VERBOSE)
        sys.exit()

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # Make output folder
    make_output_folders(data_folder)

    # Make index and hash
    make_index_and_hash(data_folder, VERBOSE=VERBOSE)

    # Map using stampy
    map_stampy(data_folder, threads=threads, VERBOSE=VERBOSE, subsample=subsample)
