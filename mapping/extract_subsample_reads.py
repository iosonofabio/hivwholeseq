#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       30/08/13
content:    Extract some reads as a subsample for fast analysis, consensus
            building, et similia.
'''
# Modules
import os
import sys
import argparse
import subprocess as sp
import numpy as np
from itertools import izip
from Bio.SeqIO.QualityIO import FastqGeneralIterator as FGI

from mapping.datasets import MiSeq_runs
from mapping.adapter_info import load_adapter_table
from mapping.filenames import get_read_filenames



# Globals
# Cluster submit
import mapping
JOBDIR = mapping.__path__[0].rstrip('/')+'/'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
JOBSCRIPT = JOBDIR+'extract_subsample_reads.py'
cluster_time = '0:59:59'
vmem = '8G'



# Functions
def fork_self(miseq_run, adaID, n_reads, VERBOSE=0, filtered=True):
    '''Fork self for each adapter ID'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+'{:02d}'.format(adaID)

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'subsam '+'{:02d}'.format(adaID),
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', miseq_run,
                 '--adaIDs', adaID,
                 '-n', n_reads,
                 '--verbose', VERBOSE,
                ]
    if not filtered:
        qsub_list.append('--raw')
    qsub_list = map(str, qsub_list)
    if VERBOSE >= 2:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def make_output_folders(data_folder, adaID, VERBOSE=0):
    '''Make the output folders if necessary'''

    subdir = data_folder+'subsample'
    if not os.path.isdir(subdir):
        os.mkdir(subdir)
        if VERBOSE:
            print 'Folder created:', subdir

    output_folder = os.path.dirname(get_read_filenames(data_folder, adaID,
                                                       subsample=True)[0])
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
        if VERBOSE:
            print 'Folder created:', output_folder


def extract_subsample(data_folder, adaID, n_reads, VERBOSE=0, filtered=True,
                      suffix=''):
    '''Extract the subsample'''
    if VERBOSE >= 1:
        print 'adapter ID: '+'{:02d}'.format(adaID)

    # Get the FASTQ files (trimmed/quality filtered if required)
    read_filenames = get_read_filenames(data_folder, adaID, subsample=False,
                                        filtered=filtered)

    # Count the number of reads
    if VERBOSE >= 2:
        print 'Counting the number of reads...',
        sys.stdout.flush()
    with open(read_filenames[0], 'r') as fh1:
        n_reads_tot = sum(1 for _ in FGI(fh1))
    if VERBOSE >= 2:
        print n_reads_tot
        sys.stdout.flush()

    # Get the random indices of the reads to store
    # (only from the middle of the file)
    ind_store = np.arange(int(0.2 * n_reads_tot), int(0.8 * n_reads_tot))
    np.random.shuffle(ind_store)
    ind_store = ind_store[:n_reads]
    ind_store.sort()

    if VERBOSE >= 2:
        print 'Random indices between '+str(ind_store[0])+' and '+str(ind_store[-1])

    # Prepare output data structures
    n_written = 0
    read1_subsample = []
    read2_subsample = []

    if VERBOSE >= 2:
        print 'Getting the reads...'
        sys.stdout.flush()

    # Iterate over the pair of files with a fast iterator
    # This is fast because we make no objects, we just save string tuples
    with open(read_filenames[0], 'r') as fh1,\
         open(read_filenames[1], 'r') as fh2:
        for i, (seq1, seq2) in enumerate(izip(FGI(fh1), FGI(fh2))):

            if VERBOSE >= 3:
                if not ((i+1) % 10000):
                    print i+1, n_written, ind_store[n_written]

            # If you hit a read, write them
            if i == ind_store[n_written]:
                read1_subsample.append(seq1)
                read2_subsample.append(seq2)
                n_written += 1

            # Break after the last one
            if n_written >= n_reads:
                break

    if VERBOSE >= 2:
        print 'done. Writing the subsample output...',
        sys.stdout.flush()

    # Write output
    out_filenames = get_read_filenames(data_folder, adaID, subsample=True,
                                       filtered=filtered, suffix=suffix)
    with open(out_filenames[0], 'w') as fh1,\
         open(out_filenames[1], 'w') as fh2:

        # Iterate over the pair of reads
        for seq1, seq2 in izip(read1_subsample, read2_subsample):
            fh1.write("@%s\n%s\n+\n%s\n" % seq1)
            fh2.write("@%s\n%s\n+\n%s\n" % seq2)

    if VERBOSE >= 2:
        print 'done.'
        print 'File written:', out_filenames[0]
        print 'File written:', out_filenames[1]
        sys.stdout.flush()



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Extract a subsample of reads')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--adaIDs', type=int, nargs='+',
                        help='Adapter ID')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--raw', action='store_true',
                        help='Use the raw reads instead of the filtered ones')
    parser.add_argument('-n', type=int, default=1000,
                        help='Subsample size')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
 
    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    VERBOSE = args.verbose
    filtered = not args.raw
    n_reads = args.n
    submit = args.submit

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # If no adapter ID is specified, iterate over all
    if not adaIDs:
        adaIDs = load_adapter_table(data_folder)['ID']

    # Iterate over adapters if needed
    for adaID in adaIDs:

        # Submit to the cluster self if requested
        if submit:
            fork_self(miseq_run, adaID, n_reads,
                      VERBOSE=VERBOSE, filtered=filtered)

        else:
            # Make output folder
            make_output_folders(data_folder, adaID, VERBOSE=VERBOSE)

            # extrace the subsample
            extract_subsample(data_folder, adaID, n_reads,
                              VERBOSE=VERBOSE, filtered=filtered)

