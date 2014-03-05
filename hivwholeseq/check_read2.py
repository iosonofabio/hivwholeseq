# vim: fdm=marker
'''
author:     Fabio Zanini
date:       25/02/14
content:    Check the demultiplexed reads.
'''
# Modules
import os
import argparse
import numpy as np
import gzip
from itertools import izip
from Bio import SeqIO

from hivwholeseq.filenames import get_read_filenames
from hivwholeseq.samples import samples
from hivwholeseq.datasets import MiSeq_runs



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Collect a few reads to perform basic checks')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28, test_tiny)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to analyze')
    parser.add_argument('--adaID', required=True,
                        help='Adapter ID to analyze (e.g. TS2)')

    args = parser.parse_args()
    seq_run = args.run
    VERBOSE = args.verbose
    maxreads = args.maxreads
    adaID = args.adaID

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # Get some reads
    fns = get_read_filenames(data_folder, adaID, gzip=True)
    with gzip.open(fns[0], 'rb') as fh1, gzip.open(fns[1], 'rb') as fh2:
        reads_iter1 = SeqIO.parse(fh1, 'fastq')
        reads_iter2 = SeqIO.parse(fh2, 'fastq')

        read_pairs = []
        inds = 20000 + np.arange(100000); np.random.shuffle(inds); inds = np.sort(inds[:10])
        ii = 0
        for irp, reads in enumerate(izip(reads_iter1, reads_iter2)):
            if irp == inds[ii]:
                read_pairs.append(reads)
                ii += 1

                if ii == len(inds):
                    break
