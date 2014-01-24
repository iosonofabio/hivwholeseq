# vim: fdm=marker
'''
author:     Fabio Zanini
date:       16/01/14
content:    Check the relation between errors, phred score, and position in the
            read. In order to call an error, the sample must be a plasmid one.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO
import pysam

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.filenames import get_mapped_phix_filename, get_phix_filename
from hivwholeseq.mapping_utils import pair_generator, convert_sam_to_bam




# Functions
def phred_errors_position(data_folder, n_cycles, maxreads=-1, posrange=[0, 1e4],
                          VERBOSE=0):
    '''Relate errors, phred score, and position in the read'''

    # DIMENSIONS: read1/read2, position in read, phred, correct/error
    counts = np.zeros((2, n_cycles / 2, 45, 2), int)

    refseq = SeqIO.read(get_phix_filename(), 'fasta')
    refmat = np.array(refseq)

    input_filename = get_mapped_phix_filename(data_folder, type='bam', filtered=True)
    if not os.path.isfile(input_filename):
        convert_sam_to_bam(input_filename)

    with pysam.Samfile(input_filename, 'rb') as bamfile:
        for irp, reads in enumerate(pair_generator(bamfile)):

            if irp == maxreads:
                if VERBOSE:
                    print 'Maximal number of read pairs reached:', maxreads
                break

            # If unmapped or unpaired, discard
            if reads[0].is_unmapped or (not reads[0].is_proper_pair) or \
               reads[1].is_unmapped or (not reads[1].is_proper_pair) or \
               (not len(reads[0].cigar)) or (not len(reads[1].cigar)):
                if VERBOSE >= 3:
                    print 'Read pair unmapped/unpaired/no CIGAR:', reads[0].qname
                continue

            for read in reads:
                count = counts[read.is_read2]
                if read.is_reverse:
                    count = count[::-1]

                pos_ref = read.pos
                pos_read = 0
                for (bt, bl) in read.cigar:
                    if bt == 1:
                        pos_read += bl
                    elif bt == 2:
                        pos_ref += bl
                    else:
                        seq = np.fromstring(read.seq[pos_read: pos_read + bl], 'S1')
                        qual = np.fromstring(read.qual[pos_read: pos_read + bl], np.int8) - 33
                        ref = refmat[pos_ref: pos_ref + bl]
                        errs = np.array(seq != ref, int)
                        poss_read = np.arange(pos_read, pos_read + bl)

                        poss_ref = np.arange(pos_ref, pos_ref + bl)
                        irg = (poss_ref >= posrange[0]) & (poss_ref < posrange[1])

                        count[[poss_read[irg], qual[irg], errs[irg]]] += 1

                        pos_read += bl
                        pos_ref += bl

    return counts


def plot_phred_errors(samplelabel, counts, ax=None, plot_ref=False,
                      color='b',
                      VERBOSE=0):
    '''Plot the error counts'''

    # 1. Plot only as a function of phred
    data = counts.sum(axis=1).sum(axis=0)
    data = 1.0 * data[:, 1] / data[:, 0]

    import matplotlib.pyplot as plt
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    ax.scatter(np.arange(len(data)), data, lw=2, label=samplelabel,
               edgecolor='none', facecolor=color)
    ax.set_xlabel('Phred score')
    ax.set_ylabel('Error rate')
    ax.set_yscale('log')
    ax.set_ylim(1e-5, 1)
    ax.set_xlim(0, 45)

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
    parser.add_argument('--runs', required=True, nargs='+',
                        help='Sequencing runs to analyze (e.g. Tue28)')
    parser.add_argument('--verbose', default=0, type=int,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to analyze')
    parser.add_argument('--posrange', type=int, nargs=2,
                        help='Range of phiX to analyze')

    args = parser.parse_args()
    seq_runs = args.runs
    VERBOSE = args.verbose
    maxreads = args.maxreads
    posrange = args.posrange

    # Use a plot for all
    from matplotlib import cm
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1)

    for i, seq_run in enumerate(seq_runs):

        dataset = MiSeq_runs[seq_run]
        data_folder = dataset['folder']
        n_cycles = dataset['n_cycles']

        counts = phred_errors_position(data_folder, n_cycles,
                                       posrange=posrange,
                                       maxreads=maxreads, VERBOSE=VERBOSE)

        plot_phred_errors(seq_run, counts, ax, plot_ref=(i == 0),
                          color=cm.jet(int(255.0 * i / len(seq_runs))),
                          VERBOSE=VERBOSE)

        plt.legend(loc='best')
