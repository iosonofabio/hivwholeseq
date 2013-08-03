#!/ebio/ag-neher/share/programs/EPD/bin/python
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
from Bio import SeqIO
from itertools import izip
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt



# Globals
VERBOSE = 1
data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/'
datafile_read1 = data_folder+'lane1_NoIndex_L001_R1_001.fastq'
datafile_adapter = data_folder+'lane1_NoIndex_L001_R2_001.fastq'
datafile_read2 = data_folder+'lane1_NoIndex_L001_R3_001.fastq'
adapters_used = {'CGATGT': 02,
                 'TGACCA': 04,
                 'CAGATC': 07,
                 'CCGTCC': 16,
                 'GTCCGC': 18,
                 'GTGAAA': 19}
adapters_table_file = 'adapters_table.dat'


# Function
def next_adapter(a1, adapters_used):
    '''Look for an adapter close to the input sequence'''
    a1 = np.array(list(a1), 'S1')
    for a in adapters_used:
        d = (a1 != np.array(list(a), 'S1')).sum()
        if d <= 1:
            return a
    return None




# Script
if __name__ == '__main__':

    # FIXME: do better with input arguments and stuff

    # Initialize adapter table file
    with open(data_folder+adapters_table_file, 'w') as f:
        f.write('\t'.join(['# adapter sequence', 'adapter ID'])+'\n')
    
    # List of found adapters
    adapters = set()

    # List of found subfolders
    g = os.walk(data_folder)
    subdirs = g.next()[1]

    # Make a default directory for unclassified reads
    if 'unclassified_reads' not in subdirs:
        os.mkdir(data_folder+'unclassified_reads')
        subdirs.append('unclassified_reads')

    # Prepare iterators
    iter_read1 = SeqIO.parse(datafile_read1, 'fastq')
    iter_adapter = SeqIO.parse(datafile_adapter, 'fastq')
    iter_read2 = SeqIO.parse(datafile_read2, 'fastq')
    iter_tot = izip(iter_read1, iter_adapter, iter_read2)

    # Iterate over all reads
    for i, (read1, adapter, read2) in enumerate(iter_tot):
        
        # Print some output
        if VERBOSE and (not ((i + 1) % 1000)):
            print i + 1

        # If the adapter is not known, add it to the list
        adapter_string = str(adapter.seq)
        adapters.add(adapter_string)

        # If the adapter does not match any know one, throw into wastebin folder
        if adapter_string in adapters_used:
            dirname = ('adapterID_'+
                       '{:02d}'.format(adapters_used[adapter_string])+
                       '/')
        else:
            dirname = 'unclassified_reads/'

        # Create a folder if not present
        if dirname.rstrip('/') not in subdirs:
            os.mkdir(data_folder+dirname)
            with open(data_folder+adapters_table_file, 'a') as f:
                f.write('\t'.join([adapter_string, str(adapters_used[adapter_string])])+'\n')
            subdirs.append(dirname.rstrip('/'))
            if VERBOSE:
                print 'Folder created:', data_folder+dirname

        # Write sequences
        with open(data_folder+dirname+'read1.fastq', 'a') as fout:
            SeqIO.write(read1, fout, 'fastq')
        with open(data_folder+dirname+'read2.fastq', 'a') as fout:
            SeqIO.write(read2, fout, 'fastq')
        if 'adapterID' not in dirname:
            with open(data_folder+dirname+'adapter.fastq', 'a') as fout:
                SeqIO.write(adapter, fout, 'fastq')

