# vim: fdm=marker
'''
author:     Fabio Zanini
date:       13/01/14
content:    Check the quality score along reads for read1 and read2.
            This rapid check is useful at the beginning, before even demultiplexing.
'''
# Modules
import os
import sys
import argparse
import numpy as np
from operator import itemgetter
from Bio import SeqIO
from itertools import izip
from Bio.SeqIO.QualityIO import FastqGeneralIterator as FGI
import pysam

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.filenames import get_demultiplex_summary_filename, get_raw_read_files, \
        get_premapped_file
from hivwholeseq.adapter_info import adapters_illumina, foldername_adapter
from hivwholeseq.fork_cluster import fork_quality_along_read as fork_self



# Functions
def quality_score_along_reads(read_len, reads_filenames,
                              maxreads=-1, VERBOSE=0):
    '''Calculate the quality score along the reads'''

    quality = [[[] for j in xrange(read_len)] for i in xrange(2)]

    # Precompute conversion table
    SANGER_SCORE_OFFSET = ord("!")
    q_mapping = dict()
    for letter in range(0, 255):
        q_mapping[chr(letter)] = letter - SANGER_SCORE_OFFSET

    # Iterate over all reads (using fast iterators)
    with open(reads_filenames[0], 'r') as fh1,\
         open(reads_filenames[1], 'r') as fh2:

        for i, reads in enumerate(izip(FGI(fh1), FGI(fh2))):

            # Stop at the maximal number of reads (for testing)
            if i == maxreads:
                if VERBOSE:
                    print 'Maximal number of read pairs reached:', maxreads
                break

            # Print some output
            if VERBOSE and (not ((i + 1) % 10000)):
                print i + 1

            for ip, read in enumerate(reads):
                for j, qletter in enumerate(read[2]):
                    quality[ip][j].append(q_mapping[qletter])

    for qual in quality:
        for qpos in qual:
            qpos.sort()

    return quality


def quality_score_along_reads_mapped(read_len, bamfilename,
                                     insertsize_range=[400, 1000],
                                     maxreads=-1, VERBOSE=0):
    '''Calculate the quality score along the reads'''
    from hivwholeseq.mapping_utils import trim_read_pair_crossoverhangs as trim_coh
    from hivwholeseq.mapping_utils import pair_generator

    quality = [[[] for j in xrange(read_len)] for i in xrange(2)]

    # Precompute conversion table
    SANGER_SCORE_OFFSET = ord("!")
    q_mapping = dict()
    for letter in range(0, 255):
        q_mapping[chr(letter)] = letter - SANGER_SCORE_OFFSET

    # Iterate over all reads (using fast iterators)
    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        for i, reads in enumerate(pair_generator(bamfile)):

            # Stop at the maximal number of reads (for testing)
            if i == maxreads:
                if VERBOSE:
                    print 'Maximal number of read pairs reached:', maxreads
                break

            # Print some output
            if VERBOSE and (not ((i + 1) % 10000)):
                print i + 1

            # Check insert size
            read = reads[reads[0].is_reverse]
            if (read.is_unmapped or (not read.is_proper_pair) or \
                (read.isize < insertsize_range[0]) or \
                (read.isize >= insertsize_range[1])):
                continue

            trim_coh(reads, trim=5, include_tests=False)

            pos_read = 0
            for read in reads:
                ip = read.is_read2
                for (bt, bl) in read.cigar:
                    if bt == 1:
                        pos_read += bl
                    elif bt == 2:
                        pass
                    elif bt == 0:
                        qualb = read.qual[pos_read: pos_read + bl]
                        poss_read = np.arange(pos_read, pos_read + bl)
                        if read.is_reverse:
                            poss_read = len(read.seq) - 1 - poss_read

                        for j, qletter in izip(poss_read, qualb):
                            quality[ip][j].append(q_mapping[qletter])

    for qual in quality:
        for qpos in qual:
            qpos.sort()

    return quality


def plot_quality_along_reads(data_folder, title, quality, VERBOSE=0, savefig=False):
    '''Plot the results of the quality scores along reads'''

    import matplotlib.pyplot as plt
    from matplotlib import cm
    fig, axs = plt.subplots(1, 2, figsize=(16, 9))
    for i, (ax, qual) in enumerate(izip(axs, quality)):
        for j, qpos in enumerate(qual):
            x = qpos
            y = np.linspace(0, 1, len(x))[::-1]
            ax.plot(x, y, color=cm.jet(int(255.0 * j / len(qual))),
                    alpha=0.5,
                    lw=2)
        ax.set_xlabel('Phred quality', fontsize=14)
        ax.set_ylabel('Fraction of bases above quality x', fontsize=14)
        ax.set_title('Read'+str(i+1), fontsize=16)
        ax.text(2, 0.03, 'blue to red: 0 to '+str(len(qual))+' base', fontsize=18)

    fig.suptitle(title, fontsize=20)

    if savefig:
        from hivwholeseq.generic_utils import mkdirs
        from hivwholeseq.filenames import get_figure_folder, \
                get_quality_along_reads_filename
        mkdirs(get_figure_folder(data_folder))
        fig.savefig(get_quality_along_reads_filename(data_folder))

    else:
        plt.ion()
        plt.show()


# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Check quality along reads')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28, test_tiny)')
    parser.add_argument('--adaID', default=None,
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to analyze')
    parser.add_argument('--submit', action='store_true', default=False,
                        help='Fork the job to the cluster via qsub')
    parser.add_argument('--no-savefig', action='store_false', dest='savefig',
                        help='Show figure instead of saving it')
    parser.add_argument('--insertsize_range', type=int, nargs=2,
                        default=(400, 1000),
                        help='Restrict to certain insert sizes')


    args = parser.parse_args()
    seq_run = args.run
    VERBOSE = args.verbose
    submit = args.submit
    maxreads = args.maxreads
    adaID = args.adaID
    savefig = args.savefig
    insertsize_range = args.insertsize_range

    # If submit, outsource to the cluster
    if submit:
        fork_self(seq_run, VERBOSE=VERBOSE, maxreads=maxreads, savefig=savefig)
        sys.exit()

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # Prepare the output data structures
    read_len = dataset['n_cycles'] / 2

    # Get quality
    if adaID is None:
        data_filenames = get_raw_read_files(dataset)
        reads_filenames = [data_filenames['read1'], data_filenames['read2']]
        quality = quality_score_along_reads(read_len, reads_filenames,
                                            maxreads=maxreads, VERBOSE=VERBOSE)

        # Plot it
        plot_quality_along_reads(data_folder, seq_run,
                                 quality, VERBOSE=VERBOSE,
                                 savefig=savefig)

    else:
        bamfilename = get_premapped_file(data_folder, adaID, type='bam')
        quality = quality_score_along_reads_mapped(read_len, bamfilename,
                                                   insertsize_range=insertsize_range,
                                                   maxreads=maxreads,
                                                   VERBOSE=VERBOSE)

        # Plot it
        plot_quality_along_reads(data_folder, seq_run+', isizes '+str(insertsize_range),
                                 quality, VERBOSE=VERBOSE,
                                 savefig=savefig)

