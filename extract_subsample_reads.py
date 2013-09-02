# vim: fdm=marker
'''
author:     Fabio Zanini
date:       30/08/13
content:    Extract some reads as a subsample for fast analysis, consensus
            building, et similia.
'''
# Modules
import sys
import argparse
import numpy as np
from itertools import izip
import Bio.SeqIO as SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator as FGI

from mapping.adapter_info import load_adapter_table
from mapping.filenames import get_read_filenames



# Globals
# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Extract a subsample of reads')
    parser.add_argument('--adaID', metavar='00', type=int, nargs='?',
                        default=0,
                        help='Adapter ID')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--raw', action='store_true',
                        help='Use the raw reads instead of the filtered ones')
    parser.add_argument('-n', type=int, default=10000,
                        help='Subsample size')
 
    args = parser.parse_args()
    adaID = args.adaID
    VERBOSE = args.verbose
    filtered = not args.raw
    n_reads = args.n

    # If no adapter ID is specified, iterate over all
    if adaID == 0:
        adaIDs = load_adapter_table(data_folder)['ID']
    else:
        adaIDs = [adaID]

    # Iterate over adapters if needed
    for adaID in adaIDs:

        if VERBOSE >= 1:
            print 'adapter ID: '+'{:02d}'.format(adaID)

        # Get the FASTQ files (trimmed/quality filtered if required)
        read_filenames = get_read_filenames(data_folder, adaID, subsample=False,
                                            filtered=filtered)

        # Get the random indices of the reads to store (between 50k and 500k only)
        ind_store = np.arange(50000, 500000)
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
            print 'Getting the reads...',
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
                                           filtered=filtered)
        with open(out_filenames[0], 'w') as fh1,\
             open(out_filenames[1], 'w') as fh2:

            # Iterate over the pair of reads
            for seq1, seq2 in izip(read1_subsample, read2_subsample):
                fh1.write("@%s\n%s\n+\n%s\n" % seq1)
                fh2.write("@%s\n%s\n+\n%s\n" % seq2)

        if VERBOSE >= 2:
            print 'done.'
            sys.stdout.flush()
