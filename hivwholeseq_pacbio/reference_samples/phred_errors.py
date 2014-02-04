# vim: fdm=marker
'''
author:     Fabio Zanini
date:       31/01/14
content:    Check the relation between errors, phred score, and position in the
            read. In order to call an error, the sample must be a plasmid one.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO
import pysam

from hivwholeseq_pacbio.datasets import data_folder_dict
from hivwholeseq_pacbio.filenames import get_premapped_file, get_reference_premap_filename



# Functions
def phred_errors_position(data_folder, samplename, n_cycles,
                          maxreads=-1, VERBOSE=0):
    '''Relate errors, phred score, and position in the read'''

    # DIMENSIONS: position in read, phred, correct/error
    counts = np.zeros((n_cycles, 81, 2), int)

    refseq = SeqIO.read(get_reference_premap_filename(data_folder, samplename), 'fasta')
    refmat = np.array(refseq)

    input_filename = get_premapped_file(data_folder, samplename, type='bam')
    with pysam.Samfile(input_filename, 'rb') as bamfile:
        for irp, read in enumerate(bamfile):

            if irp == maxreads:
                if VERBOSE:
                    print 'Maximal number of reads reached:', maxreads
                break

            if VERBOSE >= 2:
                if not ((irp+1) % 100):
                    print (irp+1)

            # If unmapped or unpaired, discard
            if read.is_unmapped or (not len(read.cigar)):
                if VERBOSE >= 3:
                    print 'Read unmapped/no CIGAR:', read.qname,
                continue

            # Exclude obvious mismaps
            if sum(bl for (bt, bl) in read.cigar if bt in (1, 2)) > 100:
                if VERBOSE >= 3:
                    print 'Obvious mismap:', read.qname
                continue


            # Collect counts
            count = counts
            if read.is_reverse:
                count = count[::-1]

            pos_ref = read.pos
            pos_read = 0
            for (bt, bl) in read.cigar:
                # Skip inserts
                if bt == 1:
                    pos_read += bl

                # Skip deletions
                elif bt == 2:
                    pos_ref += bl

                # Only mismatches are kept
                else:
                    seq = np.fromstring(read.seq[pos_read: pos_read + bl], 'S1')
                    qual = np.fromstring(read.qual[pos_read: pos_read + bl], np.int8) - 33
                    ref = refmat[pos_ref: pos_ref + bl]
                    errs = np.array(seq != ref, int)
                    
                    poss = np.arange(pos_read, pos_read + bl)
                    count[[poss, qual, errs]] += 1

                    pos_read += bl
                    pos_ref += bl

    return counts


def plot_phred_errors_position(counts, title='', label='', ax=None, plot_ref=False, VERBOSE=0):
    '''Plot the error counts'''

    # 1. Plot only as a function of phred, not of position in read
    data = counts.sum(axis=0)
    data = 1.0 * data[:, 1] / data[:, 0]

    import matplotlib.pyplot as plt
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    ax.plot(np.arange(len(data)), data, lw=2, label=label)
    ax.set_xlabel('Phred score')
    ax.set_ylabel('Error rate')
    ax.set_yscale('log')
    ax.set_ylim(1e-5, 1)
    ax.set_xlim(0, 82)
    if title:
        ax.set_title(title)

    # Reference line
    if plot_ref:
        ax.plot(np.arange(len(data)), 10**(-0.1 * np.arange(len(data))), lw=1,
                color='k', ls='--')
    
    plt.ion()
    plt.show()



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Trim and divide reads into fragments')
    parser.add_argument('--verbose', default=0, type=int,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to analyze')

    args = parser.parse_args()
    VERBOSE = args.verbose
    maxreads = args.maxreads

    # Specify the dataset
    seq_run = 'Upp23'
    samplename = 'S1'
    data_folder = data_folder_dict[seq_run]
    n_cycles = 5000

    # Use a plot for all
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1)

    counts = phred_errors_position(data_folder, samplename, n_cycles,
                                   maxreads=maxreads, VERBOSE=VERBOSE)

    plot_phred_errors_position(counts, title='PacBio NL4-3', ax=ax, plot_ref=True,
                               VERBOSE=VERBOSE)

