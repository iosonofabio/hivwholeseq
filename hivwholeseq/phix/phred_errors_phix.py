# vim: fdm=marker
'''
author:     Fabio Zanini
date:       16/01/14
content:    Check the relation between errors, phred score, and position in the
            read.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO
import pysam

from hivwholeseq.sequencing.samples import load_sequencing_run
from hivwholeseq.sequencing.filenames import get_mapped_phix_filename, get_phix_filename, \
        get_consensus_phix_filename
from hivwholeseq.utils.mapping import pair_generator, convert_sam_to_bam




# Functions
def phred_errors_position(data_folder, n_cycles, maxreads=-1, posrange=[0, 1e4],
                          VERBOSE=0, use_consensus=False):
    '''Relate errors, phred score, and position in the read'''
    if VERBOSE:
        print 'Calculating phred error by position from (filtered) mapped reads'

    # DIMENSIONS: read1/read2, position in read, phred, correct/error
    counts = np.zeros((2, n_cycles / 2, 45, 2), int)

    if use_consensus:
        ref_filename = get_consensus_phix_filename(data_folder)
    else:
        ref_filename = get_phix_filename()
    refseq = SeqIO.read(ref_filename, 'fasta')
    refmat = np.array(refseq)

    if posrange is None:
        posrange = [0, len(refseq)]

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
                      color='b', title='', VERBOSE=0):
    '''Plot the error counts'''
    if VERBOSE:
        print 'Plot phred errors'

    # Plot only as a function of phred
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
    ax.grid(True)
    if title:
        ax.set_title(title)

    # Reference line
    if plot_ref:
        ax.plot(np.arange(len(data)), 10**(-0.1 * np.arange(len(data))), lw=1,
                color='k', ls='--')
    
    plt.ion()
    plt.show()


def plot_errors_along_read(samplelabel, counts, ax=None, colors=('b', 'g'), title='',
                           logscale=True, VERBOSE=0):
    '''Plot the frequency of errors (unfiltered) along read'''
    if VERBOSE:
        print 'Plotting phix errors along read'

    # Sum all quality thresholds (that's the point)
    nu_err = counts.sum(axis=2)
    nu_err = 1.0 * nu_err[:, :, 1] / nu_err.sum(axis=2)

    for i in xrange(2):
        nu = nu_err[i]
        color = colors[i]
        label = samplelabel+', read '+str(i+1)
        ax.plot(np.arange(nu.shape[0]), nu, lw=2, color=color, label=label)
    ax.set_xlim(-1, nu.shape[0] + 1)
    ax.set_ylim(1e-4, 1.05)
    ax.set_xlabel('Position in read')
    ax.set_ylabel('Error freqs')
    ax.grid(True)
    if logscale:
        ax.set_yscale('log')
    if title:
        ax.set_title(title)

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
    parser.add_argument('--consensus', action='store_true',
                        help='Use consensus instead of reference seq')

    args = parser.parse_args()
    seq_runs = args.runs
    VERBOSE = args.verbose
    maxreads = args.maxreads
    posrange = args.posrange
    use_consensus = args.consensus

    # Use a plot for all
    from matplotlib import cm
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1)
    fig2, ax2 = plt.subplots(1, 1)


    for i, seq_run in enumerate(seq_runs):
        if VERBOSE:
            print seq_run

        dataset = load_sequencing_run(seq_run)
        data_folder = dataset['folder']
        n_cycles = dataset['cycles']

        counts = phred_errors_position(data_folder, n_cycles,
                                       posrange=posrange,
                                       maxreads=maxreads, VERBOSE=VERBOSE,
                                       use_consensus=use_consensus)

        plot_phred_errors(seq_run, counts, ax, plot_ref=(i == 0),
                          color=cm.jet(int(255.0 * i / len(seq_runs))),
                          VERBOSE=VERBOSE, title='Seq errors (PhiX) VS phred score')

        colors = [cm.jet(1.0 * (2 * i + j) / (2 * len(seq_runs))) for j in xrange(2)]
        plot_errors_along_read(seq_run, counts, ax2,
                               colors=colors,
                               VERBOSE=VERBOSE, title='Seq errors (PhiX) VS phred score')

    ax.legend(loc='upper right', fontsize=14)
    ax2.legend(loc='upper left', fontsize=14)
