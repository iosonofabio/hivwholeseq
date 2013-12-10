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
import argparse
from collections import Counter
from operator import itemgetter
from Bio import SeqIO
from itertools import izip
from Bio.SeqIO.QualityIO import FastqGeneralIterator as FGI

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.filenames import get_demultiplex_summary_filename, get_raw_read_files
from hivwholeseq.adapter_info import adapters_illumina, foldername_adapter
from hivwholeseq.fork_cluster import fork_demultiplex as fork_self



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Demultiplex HIV reads')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28, test_tiny)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--submit', action='store_true', default=False,
                        help='Fork the job to the cluster via qsub')
    parser.add_argument('--no-summary', action='store_false',
                        dest='summary',
                        help='Do not save results in a summary file')

    args = parser.parse_args()
    seq_run = args.run
    VERBOSE = args.verbose
    submit = args.submit
    summary = args.summary

    # If submit, outsource to the cluster
    if submit:
        fork_self(seq_run, VERBOSE=VERBOSE, summary=summary)
        sys.exit()

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    if summary:
        with open(get_demultiplex_summary_filename(data_folder), 'w') as f:
            f.write('Call: python demultiplex.py --run '+seq_run+' --verbose '+str(VERBOSE)+'\n')

    # Get designed adapters
    adapters_designed = tuple((adaID, adapters_illumina[adaID]) for adaID in dataset['adapters'])
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

    # Make folders for the samples
    for (adaID, s) in adapters_designed:
            dirname = foldername_adapter(adaID)
            if not os.path.isdir(data_folder+dirname):
                os.mkdir(data_folder+dirname)
                if VERBOSE:
                    print 'Folder created:', dirname

    # Make a default directory for unclassified reads
    g = os.walk(data_folder)
    subdirs = g.next()[1]
    if 'unclassified_reads' not in subdirs:
        os.mkdir(data_folder+'unclassified_reads')
        subdirs.append('unclassified_reads')

        if VERBOSE:
            print 'Folder created: unclassified reads'

    if summary:
        with open(get_demultiplex_summary_filename(data_folder), 'a') as f:
            f.write('\n')
            f.write('Folders created for samples and unclassified reads (including phix).')
            f.write('\n')

    # Get the read filenames
    data_filenames = get_raw_read_files(dataset)
    datafile_read1 = data_filenames['read1']
    datafile_read2 = data_filenames['read2']
    datafile_adapter = data_filenames['adapter']

    # Open output files
    fouts = {adaID: (open(data_folder+foldername_adapter(adaID)+'read1.fastq', 'w'),
                     open(data_folder+foldername_adapter(adaID)+'read2.fastq', 'w'))
             for adaID, _ in adapters_designed}
    fouts['unclassified'] = (open(data_folder+'unclassified_reads/read1.fastq', 'w'),
                             open(data_folder+'unclassified_reads/read2.fastq', 'w'),
                             open(data_folder+'unclassified_reads/adapter.fastq', 'w'))


    # Make sure you close the files
    try:

        # Iterate over all reads (using fast iterators)
        with open(datafile_read1, 'r') as fh1,\
             open(datafile_read2, 'r') as fh2,\
             open(datafile_adapter, 'r') as fha:

            if VERBOSE >= 3:
                print 'adaID'
                print '--------------------'

            adapters_found = Counter()
            for i, (read1, read2, adapter) in enumerate(izip(FGI(fh1), FGI(fh2),
                                                             SeqIO.parse(fha, 'fastq'))):

                # Print some output
                if VERBOSE and (not ((i + 1) % 10000)):
                    print i + 1

                # If the adapter is not known, add it to the list
                adapter_string = str(adapter.seq)
                adapters_found[adapter_string] += 1

                # If the adapter does not match any know one,
                # throw into wastebin folder
                if adapter_string not in map(itemgetter(1), adapters_designed):
                    adaID = 'unclassified'
                else:
                    adaID = dict(map(reversed, adapters_designed))[adapter_string]
            
                if VERBOSE >= 3:
                    print adaID

                # Write sequences (append to file, manual but fast)
                fouts[adaID][0].write("@%s\n%s\n+\n%s\n" % read1)
                fouts[adaID][1].write("@%s\n%s\n+\n%s\n" % read2)
                if adapter_string not in map(itemgetter(1), adapters_designed):
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
            f.write('Adapters found across all reads:\n')
            for e in adapters_found.most_common():
                f.write('\t'.join(map(str, e))+'\n')

