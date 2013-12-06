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
from Bio import SeqIO
from itertools import izip
from Bio.SeqIO.QualityIO import FastqGeneralIterator as FGI

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.filenames import get_raw_read_files
from hivwholeseq.adapter_info import adapters_LT, adapters_table_file, foldername_adapter


# Globals
# Cluster submit
import hivwholeseq
JOBDIR = hivwholeseq.__path__[0].rstrip('/')+'/'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
JOBSCRIPT = JOBDIR+'demultiplex.py'
cluster_time = '23:59:59'
vmem = '8G'



# Functions
def fork_self(miseq_run, VERBOSE=0):
    '''Fork self'''
    import subprocess as sp

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'demplx HIV',
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', miseq_run,
                 '--verbose', VERBOSE,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    sp.call(qsub_list)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Demultiplex HIV reads')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--submit', action='store_true', default=False,
                        help='Submit the job to the cluster via qsub')

    args = parser.parse_args()
    miseq_run = args.run
    VERBOSE = args.verbose
    submit = args.submit

    # If submit, outsource to the cluster
    if submit:
        fork_self(miseq_run, VERBOSE=VERBOSE)
        sys.exit()

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # Create adapter table file
    with open(data_folder+adapters_table_file, 'w') as f:
        f.write('\t'.join(['# adapter sequence', 'adapter ID'])+'\n')
    
    # List of found adapters
    adapters = set()

    # List of used adapters
    used_adapters = {adapters_LT[adaID]: adaID for adaID in dataset['adapters']}
    for adapter_string, adaID in used_adapters.iteritems():
            dirname = foldername_adapter(adaID)

            # Make subdir if necessary
            if not os.path.isdir(data_folder+dirname):
                os.mkdir(data_folder+dirname)
                
                if VERBOSE:
                    print 'Subfolder created:', dirname

            # Add entry to the adapters table
            sample = dataset['samples'][dataset['adapters'].index(adaID)]
            with open(data_folder+adapters_table_file, 'a') as f:
                f.write('\t'.join([adapter_string, str(adaID), sample])+'\n')

            if VERBOSE:
                print 'New adapter found:', adaID, adapter_string, sample

    # List of found subfolders
    g = os.walk(data_folder)
    subdirs = g.next()[1]

    # Make a default directory for unclassified reads
    if 'unclassified_reads' not in subdirs:
        os.mkdir(data_folder+'unclassified_reads')
        subdirs.append('unclassified_reads')

        if VERBOSE:
            print 'Subfolder created: unclassified reads'

    # Get the read filenames
    data_filenames = get_raw_read_files(dataset)
    datafile_read1 = data_filenames['read1']
    datafile_read2 = data_filenames['read2']
    datafile_adapter = data_filenames['adapter']

    # Open output files
    fouts = {adaID: (open(data_folder+foldername_adapter(adaID)+'read1.fastq', 'w'),
                     open(data_folder+foldername_adapter(adaID)+'read2.fastq', 'w'))
             for adaID in used_adapters.itervalues()}
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

            for i, (read1, read2, adapter) in enumerate(izip(FGI(fh1), FGI(fh2),
                                                             SeqIO.parse(fha, 'fastq'))):

                # Print some output
                if VERBOSE and (not ((i + 1) % 10000)):
                    print i + 1

                # If the adapter is not known, add it to the list
                adapter_string = str(adapter.seq)
                adapters.add(adapter_string)

                # If the adapter does not match any know one,
                # throw into wastebin folder
                if adapter_string not in used_adapters:
                    adaID = 'unclassified'
                else:
                    adaID = used_adapters[adapter_string]
            
                if VERBOSE >= 3:
                    print adaID

                # Write sequences (append to file, manual but fast)
                fouts[adaID][0].write("@%s\n%s\n+\n%s\n" % read1)
                fouts[adaID][1].write("@%s\n%s\n+\n%s\n" % read2)
                if adapter_string not in used_adapters:
                    SeqIO.write(adapter, fouts['unclassified'][2], 'fastq')

    finally:
        # Close all adaIDs
        for fout in fouts.itervalues():
            # Close both read 1 and read 2
            for fou in fout:
                fou.close()
