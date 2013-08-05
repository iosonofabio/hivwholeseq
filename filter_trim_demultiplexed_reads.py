#!/ebio/ag-neher/share/programs/EPD/bin/python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/08/13
content:    Filter and trim demultiplexed reads.
            We take either full high-quality reads or only long, high-quality
            chunks of reads.

            Note: we must take both members of a paired-end!
'''
# Modules
import os
import sys
import numpy as np
from Bio import SeqIO
from itertools import izip



# Globals
VERBOSE = 1
data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/'
adapters_table_file = 'adapters_table.dat'

# Quality threshold
phred_min = 20
block_len_min = 100


# Functions
def load_adapter_table(data_folder):
    '''Load table of adapters and samples'''
    table = np.loadtxt(data_folder+adapters_table_file,
                       dtype=[('seq', 'S6'), ('ID', int), ('sample', 'S50')],
                       ndmin=1)
    return table



# Script
if __name__ == '__main__':

    # Get a list of the good adapter IDs
    adapter_table = load_adapter_table(data_folder)

    # If an input argument is specified, do only that adapter ID, else do all
    args = sys.argv
    if len(args) > 1:
        adaIDs = [int(args[1])]
    else:
        adaIDs = adapter_table['ID']

    # Iterate over adapters (if needed)
    for adaID in adaIDs:

        # Keep log
        if VERBOSE:
            print 'adapterID: '+'{:02d}'.format(adaID)

        # Prepare temporary data structures
        read_good = np.zeros(2, bool)
        
        # Directory to read
        dirname = 'adapterID_'+'{:02d}'.format(adaID)+'/'
        
        # Scroll read files with demultiplexed raw reads
        read1_it = SeqIO.parse(data_folder+dirname+'read1.fastq', 'fastq')
        read2_it = SeqIO.parse(data_folder+dirname+'read2.fastq', 'fastq')
        reads_out = 0

        # Prepare output data structures
        with open(data_folder+dirname+'read1_filtered_trimmed.fastq', 'w') as fr1_h, \
             open(data_folder+dirname+'read2_filtered_trimmed.fastq', 'w') as fr2_h:

            # check both reads together
            for seqs12 in izip(read1_it, read2_it):
                read_good[:] = False
                reads_trimmed = []

                for ir, seq in enumerate(seqs12):
        
                    # read the phred score, make int array
                    phred = np.asarray(seq.letter_annotations['phred_quality'], int)
         
                    # get all positions above the cut-off
                    ind = np.asarray(phred >= phred_min, int)
                    # divide in blocks
                    switch = np.diff(ind).nonzero()[0] + 1
                    ind_block_start = np.insert(switch, 0, 0)
                    ind_block_end = np.append(switch, len(ind))
                    blocks_len = ind_block_end - ind_block_start
                    # get largest block
                    ind_largest_block = blocks_len.argmax()
        
                    # if there is at least one block and the largest block
                    # has at least a certain length, keep that part of the read
                    if ind.any() and (blocks_len[ind_largest_block] >= block_len_min):
                        read_good[ir] = True
                        reads_trimmed.append(seq[ind_block_start[ind_largest_block]:
                                                 ind_block_end[ind_largest_block]])
        
                # If both paired reads are good, keep them
                if read_good.all():
                    SeqIO.write(reads_trimmed[0], fr1_h, 'fastq')
                    SeqIO.write(reads_trimmed[1], fr2_h, 'fastq')
                    reads_out += 1
        
                    # Keep a log
                    if VERBOSE:
                        if not (reads_out % 10000):
                            print reads_out
