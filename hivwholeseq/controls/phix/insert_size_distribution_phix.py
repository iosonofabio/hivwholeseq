# vim: fdm=marker
'''
author:     Fabio Zanini
date:       14/10/13
content:    Quantify the distribution of insert sizes after mapping.
'''
# Modules
import os
import argparse
import pysam
import numpy as np

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.sequencing.filenames import get_mapped_phix_filename
from hivwholeseq.utils.mapping import pair_generator, convert_sam_to_bam



# Functions
def get_insert_size_distribution(data_folder, bins=None,
                                 maxreads=-1, VERBOSE=0):
    '''Get the distribution of insert sizes'''

    if maxreads > 0:
        insert_sizes = np.zeros(maxreads, np.int16)
    else:
        insert_sizes = np.zeros(1e6, np.int16)

    bamfilename = get_mapped_phix_filename(data_folder, type='bam')
    if not os.path.isfile(bamfilename):
        convert_sam_to_bam(bamfilename)

    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        # Iterate over single reads (no linkage info needed)
        n_written = 0
        for i, reads in enumerate(pair_generator(bamfile)):

            if i == maxreads:
                if VERBOSE >= 2:
                    print 'Max reads reached:', maxreads
                break
        
            # Print output
            if (VERBOSE >= 3) and (not ((i +1) % 10000)):
                print (i+1)

            # If unmapped or unpaired, mini, or insert size mini, discard
            if reads[0].is_unmapped or (not reads[0].is_proper_pair) or \
               reads[1].is_unmapped or (not reads[1].is_proper_pair):
                continue
            
            # Store insert size
            i_fwd = reads[0].is_reverse
            insert_sizes[i] = reads[i_fwd].isize
            n_written += 1

    insert_sizes = insert_sizes[:n_written]
    insert_sizes.sort()
    insert_sizes = insert_sizes[insert_sizes > 0]

    # Bin it
    if bins is None:
        h = np.histogram(insert_sizes, density=True)
    else:
        h = np.histogram(insert_sizes, bins=bins, density=True)

    return insert_sizes, h


def plot_cumulative_histogram(seq_run, insert_sizes,
                              show=False,
                              **kwargs):
    '''Plot cumulative histogram of insert sizes'''
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1)
    ax.plot(insert_sizes, np.linspace(0, 1, len(insert_sizes)), **kwargs)
    ax.set_xlabel('Insert size')
    ax.set_ylabel('Cumulative fraction')
    ax.set_title('run '+str(seq_run)+', phiX')

    plt.tight_layout()

    if show:
        plt.ion()
        plt.show()


def plot_histogram(seq_run, h,
                   show=False, savefig=False,
                   **kwargs):
    '''Plot histogram of insert sizes'''
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1)
    ax.set_title('run '+str(seq_run)+', phiX')
    x = 0.5 * (h[1][1:] + h[1][:-1])
    y = h[0]
    ax.plot(x, y, **kwargs)
    ax.set_xlabel('Insert size')
    ax.set_ylabel('Density')

    plt.tight_layout()

    if show:
        plt.ion()
        plt.show()



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Get allele counts')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--premapped', action='store_true',
                        help='Analyze premapped reads')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to analyze')
    parser.add_argument('--savefig', action='store_true',
                        help='Store figures')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    maxreads = args.maxreads
    savefig = args.savefig
    premapped = args.premapped

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # Set the bins
    bins = np.linspace(0, 1000, 100)

    # Make a single figure for the histograms
    import matplotlib.pyplot as plt
    from matplotlib import cm


    isz, h = get_insert_size_distribution(data_folder,
                                     bins=bins, maxreads=maxreads,
                                     VERBOSE=VERBOSE)
    plot_cumulative_histogram(seq_run, isz, lw=2, c='b')
    plot_histogram(seq_run, h,
                   lw=2,
                   color='b')

    plt.ion()
    plt.show()


