#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/08/13
content:    Demultiplex the FASTQ files we get into single adapter IDs (samples).
            Note: we get three files from the MiSeq: read1, read2, barcode, and
            the order of the reads is the same (we do not have to search around).
'''
# Modules
import os
import sys
import gzip
import argparse
from collections import Counter
from operator import itemgetter
from Bio import SeqIO
from itertools import izip
from Bio.SeqIO.QualityIO import FastqGeneralIterator as FGI

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.sequencing.filenames import get_demultiplex_summary_filename, get_raw_read_files, \
        get_read_filenames, get_unclassified_reads_filenames
from hivwholeseq.sequencing.adapter_info import adapters_illumina, foldername_adapter
from hivwholeseq.cluster.fork_cluster import fork_demultiplex as fork_self



# Functions
def make_output_folders(data_folder, adapters_designed, VERBOSE=0, summary=True):
    '''Make output folders for all adapters and unclassified (e.g. PhiX)'''
    from hivwholeseq.utils.generic import mkdirs

    # Make folders for the samples
    for (adaID, s) in adapters_designed:
            dirname = foldername_adapter(adaID)
            mkdirs(data_folder+dirname)
            if VERBOSE:
                print 'Folder created:', dirname

    # Make a default directory for unclassified reads
    mkdirs(data_folder+'unclassified_reads')
    if VERBOSE:
        print 'Folder created: unclassified reads'

    if summary:
        with open(get_demultiplex_summary_filename(data_folder), 'a') as f:
            f.write('\n')
            f.write('Folders created for samples and unclassified reads (including phix).')
            f.write('\n')


def get_adapters_designed(dataset, VERBOSE=0, summary=True):
    '''Get dict of the adapters designed for this experiment'''

    adapters_designed = [(adaID, adapters_illumina[adaID])
                         for adaID in dataset['adapters']]
    output = [['Designed adapters:'],
              ['adaID', 'seq', 'sample']]
    for (adaID, s) in adapters_designed:
        sample = dataset['samples'][dataset['adapters'].index(adaID)]
        output.append([adaID, s, sample])
    output = '\n'.join(map('\t'.join, output))
    if VERBOSE:
        print output
    if summary:
        with open(get_demultiplex_summary_filename(data_folder), 'a') as f:
            f.write('\n')
            f.write(output)
            f.write('\n')

    return adapters_designed


def demultiplex_reads_single_index(data_folder, data_filenames, adapters_designed,
                                   maxreads=-1, VERBOSE=0, summary=True):
    '''Demultiplex reads with single index adapters'''

    # Get the read filenames
    datafile_read1 = data_filenames['read1']
    datafile_read2 = data_filenames['read2']
    datafile_adapter = data_filenames['adapter']

    # Open output files (compressed)
    fouts = {adaID: [gzip.open(fn, 'wb', compresslevel=9)
             for fn in get_read_filenames(data_folder, adaID, gzip=True)]
             for adaID, _ in adapters_designed}

    fouts['unclassified'] = [gzip.open(fn, 'wb', compresslevel=9)
                             for fn in get_unclassified_reads_filenames(data_folder, gzip=True)]

    adapters_designed_inv = dict(map(reversed, adapters_designed))
    adapters_strings = map(itemgetter(1), adapters_designed)

    # Make sure you close the files
    try:

        # Iterate over all reads (using fast iterators)
        with gzip.open(datafile_read1, 'rb') as fh1,\
             gzip.open(datafile_read2, 'rb') as fh2,\
             gzip.open(datafile_adapter, 'rb') as fha:

            if VERBOSE >= 3:
                print 'adaID'
                print '--------------------'

            adapters_found = Counter()
            for i, (read1, read2, adapter) in enumerate(izip(FGI(fh1), FGI(fh2),
                                                             SeqIO.parse(fha, 'fastq'))):

                if i == maxreads:
                    if VERBOSE:
                        print 'Maxreads reached.'
                    break

                # Print some output
                if VERBOSE and (not ((i + 1) % 10000)):
                    print i + 1

                # If the adapter is not known, add it to the list
                adapter_string = str(adapter.seq)
                adapters_found[adapter_string] += 1

                # If the adapter does not match any know one,
                # throw into wastebin folder
                if adapter_string not in adapters_strings:
                    adaID = 'unclassified'
                else:
                    adaID = adapters_designed_inv[adapter_string]
            
                if VERBOSE >= 3:
                    print adaID

                # Write sequences (append to file, manual but fast)
                fouts[adaID][0].write("@%s\n%s\n+\n%s\n" % read1)
                fouts[adaID][1].write("@%s\n%s\n+\n%s\n" % read2)
                if adapter_string not in adapters_strings:
                    SeqIO.write(adapter, fouts['unclassified'][2], 'fastq')

    finally:
        # Close all adaIDs
        for fout in fouts.itervalues():
            # Close both read 1 and read 2 (and barcode for unclassified)
            for fou in fout:
                fou.close()

        if summary:
            with open(get_demultiplex_summary_filename(data_folder), 'a') as f:
                f.write('\n')
                f.write('Total number of reads demultiplexed: '+str(i+1)+'\n')
                f.write('Adapters found across all reads:\n')
                for e in adapters_found.most_common():
                    f.write('\t'.join(map(str, e))+'\n')


def demultiplex_reads_dual_index(data_folder, data_filenames, adapters_designed,
                                   maxreads=-1, VERBOSE=0, summary=True):
    '''Demultiplex reads with dual index adapters'''
    #FIXME: use gzipped files

    # Get the read filenames
    datafile_read1 = data_filenames['read1']
    datafile_read2 = data_filenames['read2']
    datafile_adapter1 = data_filenames['adapter1']
    datafile_adapter2 = data_filenames['adapter2']

    # Open output files
    fouts = {adaID: (open(data_folder+foldername_adapter(adaID)+'read1.fastq', 'w'),
                     open(data_folder+foldername_adapter(adaID)+'read2.fastq', 'w'))
             for adaID, _ in adapters_designed}
    fouts['unclassified'] = (open(data_folder+'unclassified_reads/read1.fastq', 'w'),
                             open(data_folder+'unclassified_reads/read2.fastq', 'w'),
                             open(data_folder+'unclassified_reads/adapter1.fastq', 'w'),
                             open(data_folder+'unclassified_reads/adapter2.fastq', 'w'),
                            )

    adapters_designed_inv = dict(map(reversed, adapters_designed))
    adapters_strings = map(itemgetter(1), adapters_designed)

    # Make sure you close the files
    try:

        # Iterate over all reads (using fast iterators)
        with open(datafile_read1, 'r') as fh1,\
             open(datafile_read2, 'r') as fh2,\
             open(datafile_adapter1, 'r') as fha1,\
             open(datafile_adapter2, 'r') as fha2:

            if VERBOSE >= 3:
                print 'adaID'
                print '--------------------'

            adapters_found = Counter()
            for i, (read1, read2,
                    adapter1,
                    adapter2) in enumerate(izip(FGI(fh1), FGI(fh2),
                                                SeqIO.parse(fha1, 'fastq'),
                                                SeqIO.parse(fha2, 'fastq'))):

                # Stop at the maximal number of reads (for testing)
                if i == maxreads:
                    if VERBOSE:
                        print 'Maximal number of read pairs reached:', maxreads
                    break

                if VERBOSE and (not ((i + 1) % 10000)):
                    print i + 1

                # If the adapter is not known, add it to the list
                adapter_string = '-'.join(map(str, [adapter1.seq, adapter2.seq]))
                adapters_found[adapter_string] += 1

                # If the adapter does not match any know one,
                # throw into wastebin folder
                if adapter_string not in adapters_strings:
                    adaID = 'unclassified'
                else:
                    adaID = adapters_designed_inv[adapter_string]
            
                if VERBOSE >= 3:
                    print adaID

                # Write sequences (append to file, manual but fast)
                fouts[adaID][0].write("@%s\n%s\n+\n%s\n" % read1)
                fouts[adaID][1].write("@%s\n%s\n+\n%s\n" % read2)
                if adapter_string not in adapters_strings:
                    SeqIO.write(adapter1, fouts['unclassified'][2], 'fastq')
                    SeqIO.write(adapter2, fouts['unclassified'][3], 'fastq')

    finally:
        # Close all adaIDs
        for fout in fouts.itervalues():
            # Close both read 1 and read 2 (and barcode for unclassified)
            for fou in fout:
                fou.close()

    if summary:
        with open(get_demultiplex_summary_filename(data_folder), 'a') as f:
            f.write('\n')
            f.write('Adapters found across all reads:\n')
            for e in adapters_found.most_common():
                f.write('\t'.join(map(str, e))+'\n')



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Demultiplex HIV reads')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28, test_tiny)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to analyze')
    parser.add_argument('--submit', action='store_true', default=False,
                        help='Fork the job to the cluster via qsub')
    parser.add_argument('--no-summary', action='store_false',
                        dest='summary',
                        help='Do not save results in a summary file')

    args = parser.parse_args()
    seq_run = args.run
    VERBOSE = args.verbose
    maxreads = args.maxreads
    submit = args.submit
    summary = args.summary

    # If submit, outsource to the cluster
    if submit:
        fork_self(seq_run, VERBOSE=VERBOSE, maxreads=maxreads, summary=summary)
        sys.exit()

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    if summary:
        with open(get_demultiplex_summary_filename(data_folder), 'w') as f:
            f.write('Call: python demultiplex.py --run '+seq_run+' --verbose '+str(VERBOSE)+'\n')

    adapters_designed = get_adapters_designed(dataset, VERBOSE=VERBOSE, summary=summary)

    make_output_folders(data_folder, adapters_designed, VERBOSE=VERBOSE,
                        summary=summary)

    data_filenames = get_raw_read_files(dataset)

    # Is it a dual index library?
    if '-' not in adapters_designed[0][0]:
        demultiplex_reads_single_index(data_folder, data_filenames, adapters_designed,
                                       maxreads=maxreads, VERBOSE=VERBOSE,
                                       summary=summary)
    else:
        demultiplex_reads_dual_index(data_folder, data_filenames, adapters_designed,
                                     maxreads=maxreads, VERBOSE=VERBOSE,
                                     summary=summary)
