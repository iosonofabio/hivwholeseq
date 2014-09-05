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
from hivwholeseq.sequencing.filenames import get_mapped_phix_filename
from hivwholeseq.sequencing.adapter_info import adapters_illumina, foldername_adapter



# Functions
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


def plot_quality_along_reads(data_folder, title, quality, VERBOSE=0):
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

    plt.ion()
    plt.show()


# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Check quality along reads')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28, test_tiny)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to analyze')
    parser.add_argument('--insertsize_range', type=int, nargs=2,
                        default=(400, 1000),
                        help='Restrict to certain insert sizes')


    args = parser.parse_args()
    seq_run = args.run
    VERBOSE = args.verbose
    maxreads = args.maxreads
    insertsize_range = args.insertsize_range

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # Prepare the output data structures
    read_len = dataset['n_cycles'] / 2

    # Get quality
    bamfilename = get_mapped_phix_filename(data_folder, type='bam')
    quality = quality_score_along_reads_mapped(read_len, bamfilename,
                                               insertsize_range=insertsize_range,
                                               maxreads=maxreads,
                                               VERBOSE=VERBOSE)

    # Plot it
    plot_quality_along_reads(data_folder, seq_run+', PhiX, isizes '+str(insertsize_range),
                             quality, VERBOSE=VERBOSE)
