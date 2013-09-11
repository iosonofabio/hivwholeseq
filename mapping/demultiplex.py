#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/08/13
content:    Demultiplex the FASTQ files we get into single adapter IDs (samples).
            Note: we get three files from the MiSeq:
                - R1 is read 2
                - R2 is adapter sequence
                - R3 is read 2
            and the order of them is the same (we do not have to search around).
'''
# Modules
import os
import sys
import argparse
from Bio import SeqIO
from itertools import izip
import numpy as np
import scipy as sp
from Bio.SeqIO.QualityIO import FastqGeneralIterator as FGI

from mapping.filenames import get_raw_read_files
from mapping.adapter_info import adapters_LT, adapters_table_file, foldername_adapter


# Globals
# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Divide HIV reads into fragments')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    VERBOSE = args.verbose

    # Create adapter table file
    with open(data_folder+adapters_table_file, 'w') as f:
        f.write('\t'.join(['# adapter sequence', 'adapter ID'])+'\n')
    
    # List of found adapters
    adapters = set()

    # List of used adapters
    used_adapters = {adapters_LT[adaID]: adaID for adaID in dataset['adapters']}
    for adaID in used_adapters.itervalues():
            dirname = foldername_adapter(adaID)

            # Make subdir if necessary
            if not os.path.isdir(data_folder+dirname):
                os.mkdir(data_folder+dirname)

            # Add entry to the adapters table
            sample = dataset['samples'][dataset['adapters'].index(adaID)]
            with open(data_folder+adapters_table_file, 'a') as f:
                f.write('\t'.join([adapter_string, str(adaID), sample])+'\n')

    # List of found subfolders
    g = os.walk(data_folder)
    subdirs = g.next()[1]

    # Make a default directory for unclassified reads
    if 'unclassified_reads' not in subdirs:
        os.mkdir(data_folder+'unclassified_reads')
        subdirs.append('unclassified_reads')

    # Get the read filenames
    data_filenames = get_raw_read_files(data_folder)
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
                if adapter_string not in adapters_used:
                    adaID = 'unclassified'
                else:
                    adaID = adapters_used[adapter_string]
            
                # Write sequences (append to file, manual but fast)
                fouts[adaID][0].write("@%s\n%s\n+\n%s\n" % read1)
                fouts[adaID][1].write("@%s\n%s\n+\n%s\n" % read2)
                if adapter_string not in adapters_used:
                    SeqIO.write(adapter, fouts[adaID][2], 'fastq')

    finally:
        for fout in fouts.itervalues():
            fout.close()
