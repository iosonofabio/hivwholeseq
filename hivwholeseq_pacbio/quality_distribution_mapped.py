# vim: fdm=marker
'''
author:     Fabio Zanini
date:       30/01/14
content:    Check the quality distribution of mapped reads (unmapped could be crap).
'''
# Modules
import os
import argparse
from itertools import izip
from Bio import SeqIO
from Bio.Seq import reverse_complement as revcom
import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#from hivwholeseq_pacbio.test_align import trim_reference, align_overlap_seqan
import seqanpy as sap
import hivwholeseq_pacbio.seqan_module.seqanpy as sap2
from hivwholeseq_pacbio.samples import samples
from hivwholeseq_pacbio.datasets import data_folder_dict
from hivwholeseq_pacbio.filenames import get_premapped_file


# Functions
def quality_histogram_single_bases(bamfilename, maxreads=-1, VERBOSE=0):
    '''Get the quality histogram of single bases'''
    qual_hist = np.zeros(100, int)
    with pysam.Samfile(bamfilename, "rb") as bamfile:
        for i, read in enumerate(bamfile):
            if i == maxreads:
                break

            if VERBOSE >= 2:
                if not ((i + 1) % 10):
                    print (i + 1)

            if read.is_unmapped:
                continue

            qual = np.fromstring(read.qual, np.int8) - 33
            for q in qual:
                qual_hist[q] += 1

    return qual_hist


def quality_average_reads(bamfilename, maxreads=-1, VERBOSE=0):
    '''Get the average quality for the whole reads'''
    quals = []
    with pysam.Samfile(bamfilename, "rb") as bamfile:
        for i, read in enumerate(bamfile):
            if i == maxreads:
                break

            if VERBOSE >= 2:
                if not ((i + 1) % 10):
                    print (i + 1)

            if read.is_unmapped:
                continue

            qual = np.fromstring(read.qual, np.int8) - 33
            quals.append(qual.mean())
    
    return quals



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Divide HIV reads into fragments')
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

    # Get quality histogram of single bases
    bamfilename = get_premapped_file(data_folder, samplename)
    qual_hist_bases = quality_histogram_single_bases(bamfilename, maxreads=maxreads,
                                                     VERBOSE=VERBOSE)

    # Get the histogram of average quality of reads
    quals_reads = quality_average_reads(bamfilename, maxreads=maxreads,
                                        VERBOSE=VERBOSE)


    # Plot single bases
    plt.figure()
    plt.plot(qual_hist_bases, lw=2, c='k')
    plt.xlabel('Phred score')
    plt.ylabel('# bases')

    plt.xlim(-1, 80)
    plt.ylim(1, maxreads * 1e3)
    plt.yscale('log')
    plt.grid()
    plt.xticks(np.arange(9) * 10)

    # Plot cumulative
    plt.figure()

    # Single bases
    plt.plot(np.arange(len(qual_hist_bases)),
             1.0 - 1.0 * np.cumsum(qual_hist_bases) / qual_hist_bases.sum(),
             c='k', lw=2, label='single bases')

    # Read average
    plt.plot(np.sort(quals_reads), 1.0 - np.linspace(0, 1, len(quals_reads)), lw=2,
             c='b', label='read average')

    plt.xlabel('Phred score')
    plt.ylabel('Fraction with q > x')
    plt.grid()

    plt.ion()
    plt.show()
