# vim: fdm=marker
'''
author:     Fabio Zanini
date:       29/06/14
content:    This script is an optional step before premapping, to trim the low-q
            end of reads in order to improve premapping (i.e. we recover more
            reads).
'''
# Modules
import sys
import os
import gzip
import time
import subprocess as sp
from itertools import izip
import argparse
import numpy as np
from Bio import SeqIO
import pysam

from hivwholeseq.samples import load_sequencing_run
from hivwholeseq.filenames import get_read_filenames, get_trim_summary_filename
from hivwholeseq.fork_cluster import fork_trim as fork_self



# Functions
def trim_reads(data_folder, adaID, VERBOSE=0, summary=True, quality=25, blocksize=10,
               minlen_read1=100, minlen_read2=50):
    '''Trim low quality at the end of reads'''
    fn_in = get_read_filenames(data_folder, adaID, gzip=True)
    fn_out = get_read_filenames(data_folder, adaID, gzip=True, trimmed=True)

    n_good = 0
    n_discarded = 0

    with gzip.open(fn_in[0], 'rb') as fin1, \
         gzip.open(fn_in[1], 'rb') as fin2, \
         gzip.open(fn_out[0], 'rb') as fout1, \
         gzip.open(fn_out[1], 'wb') as fout2:

        it1 = SeqIO.read(fin1, 'fastq')
        it2 = SeqIO.read(fin2, 'fastq')
        for irp, reads in enumerate(izip(it1, it2)):

            if VERBOSE >= 2:
                if not ((irp + 1) % 10000):
                    print irp + 1

            # Trim both reads
            trims = [trim_read(read, quality=quality, blocksize=blocksize)
                     for read in reads]

            lrs = map(len, trims)
            if (lrs[0] > minlen_read1) and (lrs[1] > minlen_read2):
                SeqIO.write(trims[0], fout1, 'fastq')
                SeqIO.write(trims[1], fout2, 'fastq')
                n_good += 1
            else:
                n_discarded += 1

    if VERBOSE:
        print 'Trim lowq ends of reads:'
        print 'Good:', n_good
        print 'Discarded:', n_discarded

    # Write summary to file
    if summary:
        with open(get_trim_summary_filename(data_folder, adaID), 'a') as f:
            f.write('\n')
            f.write('Trim low quality ends results: adaID '+adaID+'\n')
            f.write('Total:\t\t'+str(irp)+'\n')
            f.write('Good:\t\t'+str(n_good)+'\n')
            f.write('Discarded:\t'+str(n_discarded)+'\n')



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Trim loq quality end of reads')
    parser.add_argument('--run', required=True,
                        help='MiSeq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--no-summary', action='store_false', dest='summary',
                        help='Do not save results in a summary file')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    VERBOSE = args.verbose
    submit = args.submit
    summary = args.summary

    # Specify the dataset
    dataset = load_sequencing_run(seq_run)
    data_folder = dataset.folder

    # If the script is called with no adaID, iterate over all
    samples = dataset.samples
    if adaIDs is not None:
        samples = samples.loc[samples.adapter.isin(adaIDs)]
    if VERBOSE >= 2:
        print samples.index.tolist()

    # Iterate over all adaIDs
    for samplename, sample in samples.iterrows():
        adaID = str(sample.adapter)

        # Submit to the cluster self if requested
        if submit:
            fork_self(seq_run, adaID, VERBOSE=VERBOSE, threads=threads,
                      reference=refname, summary=summary)
            continue

        if summary:
            with open(get_trim_summary_filename(data_folder, adaID), 'w') as f:
                f.write('Call: python trim_reads_lowq.py --run '+seq_run+\
                        ' --adaIDs '+adaID+\
                        ' --threads '+str(threads)+\
                        ' --reference '+refname+\
                        ' --verbose '+str(VERBOSE)+'\n')

        trim_reads(data_folder, adaID, VERBOSE=VERBOSE, summary=summary)
