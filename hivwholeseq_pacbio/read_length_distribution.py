#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Richard Neher, Fabio Zanini
date:       01/08/13
content:    Plot the distribution of read lengths. On PacBio, because we do not
            have paired-end sequencing, this is almost the insert size distribution.
'''
# Modules
import sys
import os
import argparse
from Bio import SeqIO
import numpy as np
import pysam
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq_pacbio.samples import samples
from hivwholeseq_pacbio.datasets import data_folder_dict
from hivwholeseq_pacbio.filenames import get_premapped_file, \
        get_reference_premap_filename



# Functions
def get_read_length_distribution(bamfilename, maxlen=5000, maxreads=-1, VERBOSE=0):
    '''Get read length_distribution'''
    count = np.zeros(maxlen, int)

    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        for i, read in enumerate(bamfile):
            if i == maxreads:
                break

            if VERBOSE >= 2:
                if not ((i+1) % 1000):
                    print (i+1)

            if read.is_unmapped:
                if VERBOSE >= 3:
                    print 'Unmapped'
                continue

            l = len(read.seq)
            if l >= len(count):
                count_new = np.zeros(l + 1000, int)
                count_new[:len(count)] = count
                count = count_new
            count[l] += 1

    # Cut long tail
    count = count[:count.nonzero()[0][-1] + 1]

    return count


def plot_read_length_distribution(count, title):
    '''Plot the read length distribution'''

    fig, ax = plt.subplots(1, 1, figsize=(14, 6))

    # Plot cumulative
    x = np.arange(len(count))
    y = 1 - 1.0 * np.cumsum(count) / count.sum()
    
    ax.plot(x, y, lw=2, c='b')

    ax.set_xlabel('Read length [bp]')
    ax.set_ylabel('fraction of reads with length > x')
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(0, x[-1])
    ax.set_title(title)

    plt.tight_layout()
    plt.ion()
    plt.show()



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get coverage of PacBio reads')
    parser.add_argument('--run', default='Upp23',
                        help='PacBio run to analyze (e.g. Upp23)')
    parser.add_argument('--sample', required=True,
                        help='Sample to analyze (e.g. S1)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to map')

    args = parser.parse_args()
    seq_run = args.run
    samplename = args.sample
    VERBOSE = args.verbose
    maxreads = args.maxreads

    # Specify the dataset
    data_folder = data_folder_dict[seq_run]
    sample = samples.set_index('name').loc[samplename]

    # Get read length distribution
    bamfilename = get_premapped_file(data_folder, samplename)
    count = get_read_length_distribution(bamfilename, maxreads=maxreads,
                                         VERBOSE=VERBOSE)

    # Plot
    plot_read_length_distribution(count, 'PacBio read length distribution')
