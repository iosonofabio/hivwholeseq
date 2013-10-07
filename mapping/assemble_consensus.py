#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       10/09/13
content:    Assemble consensus de novo from a subsample of the reads (either
            assembling single fragments or the entire HIV genome).
'''
import os
import sys
import argparse
import subprocess as sp

from mapping.datasets import MiSeq_runs
from mapping.adapter_info import load_adapter_table
from mapping.filenames import get_read_filenames
from mapping.mapping_utils import spades_bin
from mapping.extract_subsample_reads import extract_subsample as extract_raw


# Globals
# Cluster submit
import mapping
JOBDIR = mapping.__path__[0].rstrip('/')+'/'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
JOBSCRIPT = JOBDIR+'assemble_consensus.py'
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
                 '-N', 'asb '+'{:02d}'.format(adaID),
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', miseq_run,
                 '-n', n_reads,
                 '--adaIDs', adaID,
                 '--verbose', VERBOSE,
                ]
    if not filtered:
        qsub_list.append('--raw')
    qsub_list = map(str, qsub_list)
    if VERBOSE >= 2:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def extract_subsample(data_folder, adaID, n_reads, VERBOSE=0, filtered=True):
    '''Extract a small subsample for assembly'''
    suffix = '_'+str(n_reads)
    read_filenames = get_read_filenames(data_folder, adaID, subsample=True,
                                       filtered=filtered, suffix=suffix)
    if all(map(os.path.isfile, read_filenames)):
        print >> sys.stderr, 'Subsample files found!'
    else:
        extract_raw(data_folder, adaID, n_reads,
                    VERBOSE=VERBOSE, filtered=filtered,
                    suffix=suffix)


def assemble_spades(data_folder, adaID, n_reads, VERBOSE=0, filtered=True):
    '''Assemble reads into a consensus'''
    # Input files
    suffix = '_'+str(n_reads)
    read_filenames = get_read_filenames(data_folder, adaID, subsample=True,
                                       filtered=filtered, suffix=suffix)

    # Output files
    out_folder = os.path.dirname(read_filenames[0])+'/assembly'
    call_list = [spades_bin,
                 '-t1',
                 '-m'+vmem.rstrip('G'),
                 '-k', '21,33,55,77,99,127',
                 '--careful',
                 '-o', out_folder,
                 '-1', read_filenames[0],
                 '-2', read_filenames[1]]
    if VERBOSE >= 2:
        print ' '.join(call_list)
    sp.call(call_list)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Assemble consensus de novo')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--adaIDs', type=int, nargs='+',
                        help='Adapter ID')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--raw', action='store_true',
                        help='Use the raw reads instead of the filtered ones')
    parser.add_argument('-n', type=int, default=500,
                        help='Number of reads fed to assembler')
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

        # or else, extract the sample and assemble
        else:
            extract_subsample(data_folder, adaID, n_reads,
                              VERBOSE=VERBOSE, filtered=filtered)

            assemble_spades(data_folder, adaID, n_reads,
                            VERBOSE=VERBOSE, filtered=filtered)

