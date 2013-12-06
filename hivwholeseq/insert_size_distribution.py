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
from hivwholeseq.filenames import get_mapped_filename, get_premapped_file, \
        get_insert_size_distribution_cumulative_filename, \
        get_insert_size_distribution_filename
from hivwholeseq.mapping_utils import pair_generator, convert_sam_to_bam



# Functions
def get_insert_size_distribution(data_folder, adaID, fragment, bins=None,
                                 maxreads=-1, VERBOSE=0):
    '''Get the distribution of insert sizes'''

    if maxreads > 0:
        insert_sizes = np.zeros(maxreads, np.int16)
    else:
        insert_sizes = np.zeros(1e6, np.int16)

    # Open BAM file
    if fragment == 'premapped':
        bamfilename = get_premapped_file(data_folder, adaID, type='bam')
    else:
        bamfilename = get_mapped_filename(data_folder, adaID, fragment, type='bam',
                                          filtered=True)

    # Convert from SAM if necessary
    if not os.path.isfile(bamfilename):
        convert_sam_to_bam(bamfilename)

    # Open file
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

    # Bin it
    if bins is None:
        h = np.histogram(insert_sizes, density=True)
    else:
        h = np.histogram(insert_sizes, bins=bins, density=True)

    return insert_sizes, h


def plot_cumulative_histogram(miseq_run, adaID, fragment, insert_sizes,
                              show=False, savefig=False,
                              **kwargs):
    '''Plot cumulative histogram of insert sizes'''
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1)
    ax.plot(insert_sizes, np.linspace(0, 1, len(insert_sizes)), **kwargs)
    ax.set_xlabel('Insert size')
    ax.set_ylabel('Cumulative fraction')
    ax.set_title('run '+str(miseq_run)+', adaID '+str(adaID)+', '+fragment)

    plt.tight_layout()

    if show:
        plt.ion()
        plt.show()

    if savefig:
        dataset = MiSeq_runs[miseq_run]
        data_folder = dataset['folder']
        output_filename = get_insert_size_distribution_cumulative_filename(data_folder,
                                                                           adaID,
                                                                           fragment)
        fig.savefig(output_filename)


def plot_histogram(miseq_run, adaID, fragment, h,
                   show=False, savefig=False,
                   **kwargs):
    '''Plot histogram of insert sizes'''
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1)
    x = 0.5 * (h[1][1:] + h[1][:-1])
    y = h[0]
    ax.plot(x, y, **kwargs)
    ax.set_xlabel('Insert size')
    ax.set_ylabel('Density')
    ax.set_title('run '+str(miseq_run)+', adaID '+str(adaID)+', '+fragment)

    plt.tight_layout()

    if show:
        plt.ion()
        plt.show()

    if savefig:
        dataset = MiSeq_runs[miseq_run]
        data_folder = dataset['folder']
        output_filename = get_insert_size_distribution_filename(data_folder, adaID,
                                                                fragment)
        fig.savefig(output_filename)



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Get allele counts')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to analyze')

    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    maxreads = args.maxreads

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = load_adapter_table(data_folder)['ID']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Set the bins
    bins = np.linspace(0, 1000, 100)

    # Iterate over all requested samples
    for adaID in adaIDs:
        for fragment in fragments:

            isz, h = get_insert_size_distribution(data_folder, adaID, fragment,
                                             bins=bins, maxreads=maxreads,
                                             VERBOSE=VERBOSE)
#            plot_cumulative_histogram(miseq_run, adaID, fragment, isz, lw=2, c='b')
            plot_histogram(miseq_run, adaID, fragment, h, lw=2, color='b')


