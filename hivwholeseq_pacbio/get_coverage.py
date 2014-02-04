#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/01/14
content:    Map PacBio reads to NL4-3 using stampy or custom pipeline.
'''
# Modules
import os
import argparse
from itertools import izip
from Bio import SeqIO
import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from hivwholeseq_pacbio.samples import samples
from hivwholeseq_pacbio.datasets import data_folder_dict
from hivwholeseq_pacbio.filenames import get_premapped_file, \
        get_reference_premap_filename

# Globals



# Functions
def get_coverage(bamfilename, reflen, maxreads=-1, qual_min=0, VERBOSE=0):
    '''Get the coverage from the reads'''
    cov = np.zeros(reflen, int)

    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        for i, read in enumerate(bamfile):
            if VERBOSE >= 2:
                if not ((i+1) % 1000):
                    print (i+1)
        
            if i == maxreads:
                break
        
            if read.is_unmapped:
                if VERBOSE >= 3:
                    print 'Unmapped'
                continue
        
            qual = np.fromstring(read.qual, np.int8) - 33
            pos_ref = read.pos
            pos_read = 0
            for (bt, bl) in read.cigar:

                # Matches count only above a minimal quality
                if bt == 0:
                    qualb = qual[pos_read: pos_read + bl]
                    posa = np.arange(bl)[qualb >= qual_min]
                    cov[pos_ref + posa] += 1
                    pos_ref += bl
                    pos_read += bl

                # Deletions count as covered
                elif bt == 2:
                    cov[pos_ref: pos_ref + bl] += 1
                    pos_ref += bl

    return cov


def plot_coverage(cov, title):
    '''Plot the coverage'''
    reflen = len(cov)

    fig, ax = plt.subplots(1, 1, figsize=(14, 6))
    ax.plot(cov, color='k', lw=1.5)
    ax.set_xlabel('Position [bp]')
    ax.set_ylabel('Coverage')

    ax.set_xlim(-50, reflen + 50) 
    ax.set_ylim(-10, ax.get_ylim()[1])
    ax.set_title(title)

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
    parser.add_argument('--qual_min', type=int, default=0,
                        help='Minimal Phred quality score to count')

    args = parser.parse_args()
    seq_run = args.run
    samplename = args.sample
    VERBOSE = args.verbose
    maxreads = args.maxreads
    qual_min = args.qual_min

    # Specify the dataset
    data_folder = data_folder_dict[seq_run]
    sample = samples.set_index('name').loc[samplename]

    # Get NL4-3 reference
    refseq = SeqIO.read(get_reference_premap_filename(data_folder, samplename), 'fasta')
    refm = np.array(refseq)

    # Get coverage
    bamfilename = get_premapped_file(data_folder, samplename)
    cov = get_coverage(bamfilename, len(refseq), maxreads=maxreads,
                       VERBOSE=VERBOSE,
                       qual_min=qual_min)

    # Plot
    plot_coverage(cov, 'PacBio coverage: '+samplename)
