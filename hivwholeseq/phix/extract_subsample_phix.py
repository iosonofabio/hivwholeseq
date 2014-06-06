# vim: fdm=marker
'''
author:     Fabio Zanini
date:       30/05/14
content:    Extract a few reds to blast and similia.
'''
# Modules
import os
import sys
import argparse
import gzip
import numpy as np
from Bio import SeqIO
from itertools import izip
from Bio.SeqIO.QualityIO import FastqGeneralIterator as FGI

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.filenames import get_phix_filename, \
        get_unclassified_reads_filenames



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Map reads to PhiX')
    parser.add_argument('--run', required=True,
                        help='MiSeq run to analyze (e.g. Tue28)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to analyze')

    args = parser.parse_args()
    seq_run = args.run
    VERBOSE = args.verbose
    maxreads = args.maxreads

    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    reads_filenames = get_unclassified_reads_filenames(data_folder, gzip=True)[:2]
    out_filenames = [r.replace('.fastq', '_subsample.fastq') for r in reads_filenames]
    n_reads_filename = '/'.join(reads_filenames[0].split('/')[:-1]+['n_reads.dat'])

    if reads_filenames[0][-3:] == '.gz':
        openf = gzip.open
        file_readmode = 'rb'
        file_writemode = 'wb'
    else:
        openf = open
        file_readmode = 'r'
        file_writemode = 'w'

    # Iterate over all reads (using fast iterators)
    with openf(reads_filenames[0], file_readmode) as fh1, \
         openf(reads_filenames[1], file_readmode) as fh2, \
         openf(out_filenames[0], file_writemode) as fo1, \
         openf(out_filenames[1], file_writemode) as fo2:
    
        if VERBOSE:
            print 'Getting number of reads',
        if not os.path.isfile(n_reads_filename):
            n_reads = sum(1 for read in FGI(fh1))
            fh1.rewind()
            with open(n_reads_filename, 'w') as fnr:
                fnr.write(str(n_reads)+'\n')
        else:
            with open(n_reads_filename, 'r') as fnr:
                n_reads = int(fnr.readline().rstrip('\n'))

        if VERBOSE:
            print n_reads

        inds = np.arange(n_reads)
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

            fo1.write("@%s\n%s\n+\n%s\n" % reads[0])
            fo2.write("@%s\n%s\n+\n%s\n" % reads[1])

            indi += 1
            if indi == maxreads:
                if VERBOSE:
                    print 'Maximal number of read pairs reached:', maxreads
                break


