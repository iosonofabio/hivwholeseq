# vim: fdm=marker
'''
author:     Fabio Zanini
date:       16/09/13
content:    Study the read length distribution at the end of the mapping pipeline.
'''
# Modules
import os
import argparse
import pysam
import numpy as np
import matplotlib.cm as cm

# Matplotlib parameters
import matplotlib
params = {'axes.labelsize': 20, 
          'text.fontsize': 20,
          'legend.fontsize': 18,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16,
          'text.usetex': False}
matplotlib.rcParams.update(params)
import matplotlib.pyplot as plt

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.sequencing.adapter_info import load_adapter_table
from hivwholeseq.miseq import read_types
from hivwholeseq.sequencing.filenames import get_mapped_filename
from hivwholeseq.mapping_utils import convert_sam_to_bam



def get_read_lengths(data_folder, adaID, fragment, VERBOSE=0, maxreads=-1):
    '''Get the read lengths'''

    # Lengths from 1 to 250
    lengths = np.zeros((len(read_types), 250), int)

    # Open BAM file
    # Note: the reads should already be filtered of unmapped stuff at this point
    bamfilename = get_mapped_filename(data_folder, adaID, fragment, type='bam',
                                      filtered=True)
    if not os.path.isfile(bamfilename):
        convert_sam_to_bam(bamfilename)
    with pysam.Samfile(bamfilename, 'rb') as bamfile:

        # Iterate over single reads (no linkage info needed)
        for i, read in enumerate(bamfile):

            # Max number of reads
            if i == maxreads:
                if VERBOSE >= 2:
                    print 'Max reads reached:', maxreads
                break
        
            # Print output
            if (VERBOSE >= 3) and (not ((i +1) % 10000)):
                print (i+1)
        
            # Divide by read 1/2 and forward/reverse
            js = 2 * read.is_read2 + read.is_reverse

            # Increment counter
            lengths[js, read.rlen - 1] += 1

            # Note: we do not delve into CIGARs because the reads are trimmed

    return lengths


def plot_read_length_distribution(adaID, fragment, length):
    '''Plot the distribution'''
    fig, ax = plt.subplots(1, 1)
    for irt, read_type in enumerate(read_types):
        color = cm.jet(int(255.0 * irt / len(read_types)))
        #ax.plot(np.arange(1, 251), length[irt], lw=1.5, c=color)
        ax.scatter(np.arange(1, 251), length[irt], s=50, c=color, label=read_type)
    ax.set_xlabel('Read length [bases]')
    ax.set_ylabel('# reads')
    ax.set_title('adaID '+adaID+', '+fragment)
    ax.set_yscale('log')
    ax.legend(loc='lower right', fontsize=10)
    ax.set_ylim(ymin=1)
    ax.set_xlim(xmin=-5)
    plt.tight_layout()

    ax2 = fig.add_axes([0.2, 0.58, 0.27, 0.27])
    for irt, read_type in enumerate(read_types):
        color = cm.jet(int(255.0 * irt / len(read_types)))
        ax2.scatter(np.arange(1, 251), length[irt], s=50, c=color, label=read_type)
    ax2.set_ylim(ymin=0)
    ax2.yaxis.get_major_formatter().set_powerlimits((0,1))
    ax2.set_xticks((0, 250))
    plt.show()


def plot_read_length_distribution_cumulative(adaID, fragment, length):
    '''Plot the cumulative distribution'''
    fig, ax = plt.subplots(1, 1)
    labss = {'read1 f': 'read1 fwd', 'read1 r': 'read1 rev',
             'read2 f': 'read2 fwd', 'read2 r': 'read2 rev'}
    for irt, read_type in enumerate(read_types):
        color = cm.jet(int(255.0 * irt / len(read_types)))
        ax.plot(np.arange(1, 251), 1.0 - 1.0 * np.cumsum(length[irt]) / length[irt].sum(),
                lw=1.5, c=color,
                label=labss[read_type])
    ax.set_xlabel('Read length [bases]')
    ax.set_ylabel('fraction of reads longer than x')
    #ax.set_title('adaID '+adaID+', '+fragment)
    ax.legend(loc='lower left', fontsize=18)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(-5, 255)
    plt.tight_layout()
    plt.show()






# Script
if __name__ == '__main__':


    # Input arguments
    parser = argparse.ArgumentParser(description='Study minor allele frequency')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('-n', type=int, default=-1,
                        help='Number of reads analyzed')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    n_reads = args.n
    VERBOSE = args.verbose

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
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

    # Iterate over all requested samples
    lengths_all = {}
    for adaID in adaIDs:
        for fragment in fragments:
            lengths_all[(adaID, fragment)] = length = \
                    get_read_lengths(data_folder, adaID, fragment, VERBOSE=VERBOSE,
                                     maxreads=n_reads)

            # Plot it
            #plot_read_length_distribution(adaID, fragment, length)
            plot_read_length_distribution_cumulative(adaID, fragment, length)

