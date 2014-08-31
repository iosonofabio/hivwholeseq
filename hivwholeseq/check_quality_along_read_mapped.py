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
import gzip
import numpy as np
from operator import itemgetter
from Bio import SeqIO
from itertools import izip
from Bio.SeqIO.QualityIO import FastqGeneralIterator as FGI
import pysam

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.filenames import get_demultiplex_summary_filename, get_raw_read_files, \
        get_premapped_filename, get_read_filenames
from hivwholeseq.adapter_info import adapters_illumina, foldername_adapter
from hivwholeseq.mapping_utils import extract_mapped_reads_subsample_open



# Functions
def quality_score_along_reads(read_len, reads_filenames,
                              skipreads=0,
                              randomreads=False,
                              maxreads=-1, VERBOSE=0):
 != 'nan'    '''Calculate the quality score along the reads'''

    quality = [[[] for j in xrange(read_len)] for i in xrange(2)]

    # Precompute conversion table
    SANGER_SCORE_OFFSET = ord("!")
    q_mapping = dict()
    for letter in range(0, 255):
        q_mapping[chr(letter)] = letter - SANGER_SCORE_OFFSET

    if reads_filenames[0][-3:] == '.gz':
        openf = gzip.open
        file_readmode = 'rb'
    else:
        openf = open
        file_readmode = 'r'

    # Iterate over all reads (using fast iterators)
    with openf(reads_filenames[0], file_readmode) as fh1, \
         openf(reads_filenames[1], file_readmode) as fh2:
    
        if randomreads:
            if VERBOSE:
                print 'Getting number of reads',
            n_reads = sum(1 for read in FGI(fh1))
            fh1.rewind()
            if VERBOSE:
                print n_reads

            inds = np.arange(skipreads, n_reads)
            np.random.shuffle(inds)
            inds = inds[:maxreads]
            inds.sort()
            indi = 0
            if VERBOSE:
                print 'Random indices from ', inds[0], 'to', inds[-1]

            for (i, reads) in enumerate(izip(FGI(fh1), FGI(fh2))):
                if VERBOSE and (not ((i + 1) % 10000)):
                    print i + 1

                if (i != inds[indi]):
                    continue

                for ip, read in enumerate(reads):
                    for j, qletter in enumerate(read[2]):
                        quality[ip][j].append(q_mapping[qletter])

                indi += 1
                if indi == maxreads:
                    if VERBOSE:
                        print 'Maximal number of read pairs reached:', maxreads
                    break



        else:

            for (i, reads) in enumerate(izip(FGI(fh1), FGI(fh2))):
                if i < skipreads:
                    continue

                if i == skipreads + maxreads:
                    if VERBOSE:
                        print 'Maximal number of read pairs reached:', maxreads
                    break

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
                                     skipreads=0,
                                     maxreads=-1,
                                     randomreads=True,
                                     VERBOSE=0):
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
        if not randomreads:
            reads_all = []
            for i, read_pair in enumerate(pair_generator(bamfile)):
                if i < skipreads:
                    continue
    
                if i == skipreads + maxreads:
                    if VERBOSE:
                        print 'Maximal number of read pairs reached:', maxreads
                    break
    
                if VERBOSE and (not ((i + 1) % 10000)):
                    print i + 1

                reads_all.append(read_pair)

        else:
            reads_all = extract_mapped_reads_subsample_open(bamfile, maxreads,
                                                            VERBOSE=VERBOSE)

        print len(reads_all)

        for reads in reads_all:

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


def plot_quality_along_reads(data_folder, adaID, title, quality, VERBOSE=0, savefig=False):
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
        fig_folder = get_figure_folder(data_folder, adaID)
        fig_filename = get_quality_along_reads_filename(data_folder, adaID)
        mkdirs(fig_folder)
        fig.savefig(fig_filename)

    else:
        plt.tight_layout()
        plt.ion()
        plt.show()


def plot_cuts_quality_along_reads(data_folder, adaID, title, quality, VERBOSE=0, savefig=False):
    '''Plot some cuts of the quality along the read'''
    from scipy.stats import percentileofscore as pof
    import matplotlib.pyplot as plt
    from matplotlib import cm
    fig, axs = plt.subplots(1, 2, figsize=(14, 8))
    qthreshs = [10, 20, 30, 35]
    for i, (ax, qual) in enumerate(izip(axs, quality)):
        for j, qthresh in enumerate(qthreshs):
            x = np.arange(len(qual))
            y = np.array([100 - pof(qual[k], qthresh) for k in xrange(len(qual))])
            ax.plot(x, y, color=cm.jet(int(255.0 * j / len(qthreshs))),
                    alpha=0.8,
                    lw=2,
                    label='Q = '+str(qthresh))
        ax.set_xlabel('Position [bp]', fontsize=14)
        ax.set_ylabel('Percentage of bases above quality x', fontsize=14)
        ax.set_title('Read'+str(i+1), fontsize=16)
        ax.set_ylim(-1, 101)
        ax.set_xlim(-1, len(qual) + 1)
        ax.legend(loc='best')

    fig.suptitle(title, fontsize=20)

    if savefig:
        from hivwholeseq.generic_utils import mkdirs
        from hivwholeseq.filenames import get_figure_folder, \
                get_quality_along_reads_filename
        fig_folder = get_figure_folder(data_folder, adaID)
        fig_filename = get_quality_along_reads_filename(data_folder, adaID, simple=True)
        mkdirs(fig_folder)
        fig.savefig(fig_filename)

    else:
        plt.tight_layout()
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
    parser.add_argument('--skipreads', type=int, default=0,
                        help='Skip the first reads of the file')    
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to analyze')
    parser.add_argument('--submit', action='store_true', default=False,
                        help='Fork the job to the cluster via qsub')
    parser.add_argument('--no-savefig', action='store_false', dest='savefig',
                        help='Show figure instead of saving it')
    parser.add_argument('--mapped', action='store_true',
                        help='Analyze mapped reads')
    parser.add_argument('--random', action='store_true',
                        help='Take random reads instead of consecutive ones')
    parser.add_argument('--insertsize_range', type=int, nargs=2,
                        default=(400, 1000),
                        help='Restrict to certain insert sizes')
    parser.add_argument('--plotfull', action='store_true',
                        help='Make full plots')

    args = parser.parse_args()
    seq_run = args.run
    VERBOSE = args.verbose
    submit = args.submit
    skipreads = args.skipreads
    maxreads = args.maxreads
    adaID = args.adaID
    savefig = args.savefig
    use_mapped = args.mapped
    use_random = args.random
    insertsize_range = args.insertsize_range
    plotfull = args.plotfull

    if submit:
        raise ValueError('Not implemented!')
        #fork_self(seq_run, VERBOSE=VERBOSE, maxreads=maxreads, savefig=savefig)
        #sys.exit()

    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']
    read_len = dataset['n_cycles'] / 2

    if (adaID is None) and use_mapped:
        raise ValueError('Choose an adaID if you want to study mapped reads')

    if not use_mapped:
        if adaID is None:
            data_filenames = get_raw_read_files(dataset)
            reads_filenames = [data_filenames['read1'], data_filenames['read2']]
            title = seq_run
        else:
            reads_filenames = get_read_filenames(data_folder, adaID, gzip=True)
            if not os.path.isfile(reads_filenames[0]):
                reads_filenames = get_read_filenames(data_folder, adaID, gzip=False)
            title = seq_run+', '+adaID


        quality = quality_score_along_reads(read_len, reads_filenames,
                                            skipreads=skipreads,
                                            randomreads=use_random,
                                            maxreads=maxreads, VERBOSE=VERBOSE)

    else:
        bamfilename = get_premapped_filename(data_folder, adaID, type='bam')
        title = seq_run+', isizes '+str(insertsize_range)
        quality = quality_score_along_reads_mapped(read_len, bamfilename,
                                                   insertsize_range=insertsize_range,
                                                   skipreads=skipreads,
                                                   maxreads=maxreads,
                                                   randomreads=use_random,
                                                   VERBOSE=VERBOSE)

    plot_cuts_quality_along_reads(data_folder, adaID, title,
                                  quality, VERBOSE=VERBOSE,
                                  savefig=savefig)

    if plotfull:
        plot_quality_along_reads(data_folder, adaID, title,
                                 quality, VERBOSE=VERBOSE,
                                 savefig=savefig)

