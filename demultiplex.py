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
import matplotlib.pyplot as plt

from mapping.filenames import get_raw_read_files
from adapter_info import adapters_LT, adapters_table_file


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

    # Iterate over all reads
    iter_tot = izip(SeqIO.parse(datafile_read1, 'fastq'),
                    SeqIO.parse(datafile_read2, 'fastq'),
                    SeqIO.parse(datafile_adapter, 'fastq'))
    for i, (read1, read2, adapter) in enumerate(iter_tot):
        
        # Print some output
        if VERBOSE and (not ((i + 1) % 1000)):
            print i + 1

        # If the adapter is not known, add it to the list
        adapter_string = str(adapter.seq)
        adapters.add(adapter_string)

        # If the adapter does not match any know one, throw into wastebin folder
        if adapter_string in adapters_used:
            dirname = 'adapterID_'+'{:02d}'.format(adapters_used[adapter_string])+'/'
        else:
            dirname = 'unclassified_reads/'

        # Create a subsample folder and split there the first reads for consensus building
        if i < reads_subsample:
            data_folders = [data_folder, data_folder_subsample]
        elif i == reads_subsample:
            print 'Subsample done!'
            sys.stdout.flush()
        else:
            data_folders = [data_folder]

        # Create a folder for the adapter if not present
        if dirname.rstrip('/') not in subdirs[data_folder]:
            os.mkdir(data_folder+dirname)
            with open(data_folder+adapters_table_file, 'a') as f:
                f.write('\t'.join([adapter_string, str(adapters_used[adapter_string])])+'\n')
            subdirs[data_folder].append(dirname.rstrip('/'))
            if VERBOSE:
                print 'Folder created:', data_folder+dirname
    
        # Write sequences (append to file)
        with open(data_folder+dirname+'read1.fastq', 'a') as fout:
            SeqIO.write(read1, fout, 'fastq')
        with open(data_folder+dirname+'read2.fastq', 'a') as fout:
            SeqIO.write(read2, fout, 'fastq')
        if 'adapterID' not in dirname:
            with open(data_folder+dirname+'adapter.fastq', 'a') as fout:
                SeqIO.write(adapter, fout, 'fastq')
    

    # Note: the adapter table must be enriched of descriptions by hand afterwards
