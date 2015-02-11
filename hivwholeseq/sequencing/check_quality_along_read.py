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

from hivwholeseq.sequencing.filenames import get_read_filenames
from hivwholeseq.cluster.fork_cluster import fork_quality_along_read as fork_self
from hivwholeseq.utils.mapping import extract_mapped_reads_subsample_open
from hivwholeseq.sequencing.samples import load_sequencing_run



# Functions
def quality_score_along_reads(read_len, reads_filenames,
                              skipreads=0,
                              randomreads=False,
                              maxreads=-1, VERBOSE=0):
    '''Calculate the quality score along the reads'''

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
                sys.stdout.flush()
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
        from hivwholeseq.utils.generic import mkdirs
        from hivwholeseq.sequencing.filenames import get_figure_folder, \
                get_quality_along_reads_filename
        fig_folder = get_figure_folder(data_folder, adaID)
        fig_filename = get_quality_along_reads_filename(data_folder, adaID)
        mkdirs(fig_folder)
        fig.savefig(fig_filename)

    else:
        plt.tight_layout()
        plt.ion()
        plt.show()


def plot_cuts_quality_along_reads(data_folder, adaID, quality, title='',
                                  VERBOSE=0, savefig=False):
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

    if title:
        fig.suptitle(title, fontsize=20)

    if savefig:
        from hivwholeseq.utils.generic import mkdirs
        if savefig == True:
            from hivwholeseq.sequencing.filenames import get_figure_folder, \
                    get_quality_along_reads_filename
            fig_folder = get_figure_folder(data_folder, adaID)
            fig_filename = get_quality_along_reads_filename(data_folder, adaID, simple=True)
        elif isinstance(savefig, basestring):
            import os
            fig_filename = savefig
            fig_folder = os.path.dirname(fig_filename)

        else:
            raise ValueError('savefig must be a bool or a figure filename (string)')

        mkdirs(fig_folder)
        fig.savefig(fig_filename)

    else:
        plt.tight_layout()
        plt.ion()
        plt.show()



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Check quality along reads',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28, test_tiny)')
    parser.add_argument('--adaID', required=True,
                        help='Adapter ID to analyze (e.g. TS2)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to analyze')
    parser.add_argument('--submit', action='store_true', default=False,
                        help='Fork the job to the cluster via qsub')
    parser.add_argument('--no-savefig', action='store_false', dest='savefig',
                        help='Show figure instead of saving it')

    args = parser.parse_args()
    seq_run = args.run
    VERBOSE = args.verbose
    submit = args.submit
    maxreads = args.maxreads
    adaID = args.adaID
    savefig = args.savefig

    if submit:
        fork_self(seq_run, VERBOSE=VERBOSE, maxreads=maxreads, savefig=savefig)
        sys.exit()

    dataset = load_sequencing_run(seq_run)
    data_folder = dataset.folder
    read_len = dataset.cycles // 2

    reads_filenames = get_read_filenames(data_folder, adaID, gzip=True)
    if not os.path.isfile(reads_filenames[0]):
        reads_filenames = get_read_filenames(data_folder, adaID, gzip=False)
    title = seq_run+', '+adaID

    quality = quality_score_along_reads(read_len, reads_filenames,
                                        randomreads=(maxreads >= 1),
                                        maxreads=maxreads, VERBOSE=VERBOSE)

    plot_cuts_quality_along_reads(data_folder, adaID,
                                  quality,
                                  title=title,
                                  VERBOSE=VERBOSE,
                                  savefig=savefig)

    #if plotfull:
    #    plot_quality_along_reads(data_folder, adaID, title,
    #                             quality, VERBOSE=VERBOSE,
    #                             savefig=savefig)

