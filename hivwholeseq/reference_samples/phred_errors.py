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
from hivwholeseq.sequencing.filenames import get_premapped_filename, get_reference_premap_filename, \
        get_consensus_filename, get_mapped_filename
from hivwholeseq.utils.mapping import pair_generator, convert_sam_to_bam




# Functions
def phred_errors_position(data_folder, adaID, n_cycles,
                          fragment='premapped',
                          maxreads=-1, VERBOSE=0):
    '''Relate errors, phred score, and position in the read'''

    # DIMENSIONS: read1/read2, position in read, phred, correct/error
    counts = np.zeros((2, n_cycles / 2, 45, 2), int)

    if fragment == 'premapped':
        refseq = SeqIO.read(get_reference_premap_filename(data_folder, adaID), 'fasta')
        input_filename = get_premapped_filename(data_folder, adaID, type='bam')
    else:
        refseq = SeqIO.read(get_consensus_filename(data_folder, adaID, fragment), 'fasta')
        input_filename = get_mapped_filename(data_folder, adaID, fragment, type='bam',
                                             filtered=False)

    refmat = np.array(refseq)

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
                        
                        poss = np.arange(pos_read, pos_read + bl)
                        count[[poss, qual, errs]] += 1

                        pos_read += bl
                        pos_ref += bl

    return counts


def plot_phred_errors_position(samplelabel, counts, ax=None, plot_ref=False, VERBOSE=0):
    '''Plot the error counts'''

    # 1. Plot only as a function of phred
    data = counts.sum(axis=1).sum(axis=0)
    data = 1.0 * data[:, 1] / data[:, 0]

    import matplotlib.pyplot as plt
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    ax.plot(np.arange(len(data)), data, lw=2, label=samplelabel)
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
    parser.add_argument('--run', required=True,
                        help='Sequencing run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--fragment', default='premapped',
                        help='What fragment to use/premapped')
    parser.add_argument('--verbose', default=0, type=int,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to analyze')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    fragment = args.fragment
    VERBOSE = args.verbose
    maxreads = args.maxreads

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # Set the number of cycles of the kit (for trimming adapters in short inserts)
    n_cycles = dataset['n_cycles']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = MiSeq_runs[seq_run]['adapters']


    # Use a plot for all
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1)

    # Iterate over all adaIDs
    for i, adaID in enumerate(adaIDs):

        samplename = dataset['samples'][dataset['adapters'].index(adaID)]

        counts = phred_errors_position(data_folder, adaID, n_cycles,
                                       fragment=fragment,
                                       maxreads=maxreads, VERBOSE=VERBOSE)

        plot_phred_errors_position(samplename, counts, ax, plot_ref=(i == 0),
                                   VERBOSE=VERBOSE)

    plt.legend(loc='best')
