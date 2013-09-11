#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/08/13
content:    Filter and trim demultiplexed reads.
            We take either full high-quality reads or only long, high-quality
            chunks of reads.

            Note: we must take both members of a paired-end!
'''
# Modules
import os
import sys
import argparse
import numpy as np
from Bio import SeqIO
from itertools import izip

from mapping.adapter_info import load_adapter_table, foldername_adapter
from mapping.filenames import get_read_filenames, get_read_unpaired_filename



# Globals
# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']

# Quality threshold
phred_min = 20
block_len_min = 100


# Cluster submit
import mapping
JOBDIR = mapping.__path__[0].rstrip('/')+'/'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
JOBSCRIPT = JOBDIR+'filter_trim_demultiplexed_reads.py'
cluster_time = '0:59:59'
vmem = '8G'



# Functions
def fork_self(data_folder, adaID, fragment, VERBOSE=0):
    '''Fork self for each adapter ID'''
    import subprocess as sp

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'exm_'+'{:02d}'.format(adaID),
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--adaIDs', adaID,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    sp.call(qsub_list)
    

def filter_trim_reads(data_folder, adaID, VERBOSE=0, subsample=False):
    '''Filter the reads to good chunks, based on phred quality'''

    # Prepare temporary data structures
    read_good = np.zeros(2, bool)
    
    # Directory to read
    dirname = foldername_adapter(adaID)
    
    # Scroll read files with demultiplexed raw reads
    readfiles = get_read_filenames(data_folder, adaID, subsample=subsample,
                                   filtered=False)
    read1_it = SeqIO.parse(readfiles[0], 'fastq')
    read2_it = SeqIO.parse(readfiles[1], 'fastq')

    # Open output files
    outfiles = get_read_filenames(data_folder, adaID, subsample=subsample,
                                  filtered=True)
    outfile_unpaired = get_read_unpaired_filename(data_folder, adaID,
                                                  subsample=subsample)
    reads_out = 0
    read_unpaired = 0
    reads_missed = 0
    with open(outfiles[0], 'w') as fr1_h, \
         open(outfiles[1], 'w') as fr2_h, \
         open(outfile_unpaired, 'w') as fru_h:

        # Iterate over read pairs
        for reads in izip(read1_it, read2_it):

            read_good[:] = False
            reads_trimmed = []

            for ir, seq in enumerate(reads):
    
                # read the phred score, make int array
                phred = np.asarray(seq.letter_annotations['phred_quality'], int)
     
                # get all positions above the cut-off
                ind = np.asarray(phred >= phred_min, int)
                # divide in blocks
                switch = np.diff(ind).nonzero()[0] + 1
                ind_block_start = np.insert(switch, 0, 0)
                ind_block_end = np.append(switch, len(ind))
                blocks_len = ind_block_end - ind_block_start
                # get largest block
                ind_largest_block = blocks_len.argmax()
    
                # if there is at least one block and the largest block
                # has at least a certain length, keep that part of the read
                if ind.any() and (blocks_len[ind_largest_block] >= block_len_min):
                    read_good[ir] = True
                    reads_trimmed.append(seq[ind_block_start[ind_largest_block]:
                                             ind_block_end[ind_largest_block]])
    
            # If both paired reads are good, keep them
            if read_good.all():
                SeqIO.write(reads_trimmed[0], fr1_h, 'fastq')
                SeqIO.write(reads_trimmed[1], fr2_h, 'fastq')
                reads_out += 1

            # If only one is good, the read is essentially not paired, store
            # it somewhere
            elif read_good.sum() == 1:
                SeqIO.write(reads_trimmed[0], fru_h, 'fastq')
                SeqIO.write(reads_trimmed[1], fru_h, 'fastq')
                reads_unpaired += 1

                # Keep a log
                if VERBOSE >= 3:
                    if not (reads_out % 10000):
                        print reads_out, reads_unpaired

            else:
                reads_missed += 1


    # Keep a log
    if VERBOSE >= 2:
        print 'Read pairs good, unpaired, missed:',
        print reads_out, reads_unpaired, reads_missed





# Script
if __name__ == '__main__':

    # Parse input arguments
    parser = argparse.ArgumentParser(description='Filter & trim demultiplexed reads.')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--subsample', action='store_true',
                        help='Apply only to a subsample of the reads')
    parser.add_argument('--verbose', type=int, default=0,
                        help=('Verbosity level [0-3]'))
    parser.add_argument('--submit', action='store_true', default=False,
                        help='Submit the job to the cluster via qsub')

    args = parser.parse_args()
    adaIDs = list(args.adaIDs)
    VERBOSE = args.verbose
    submit = args.submit

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = load_adapter_table(data_folder)['ID']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # Iterate over adapters
    for adaID in adaIDs:

        # Submit to the cluster self if requested
        if submit:
            fork_self(data_folder, adaID, VERBOSE=VERBOSE, subsample=subsample)
            continue

        # Filter ad trim
        filter_trim_reads(data_folder, adaID, VERBOSE=VERBOSE, subsample=subsample)
