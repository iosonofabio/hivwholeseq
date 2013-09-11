# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/08/13
content:    Take demultiplexed reads and pick a subsample to build a good consensus.

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
from map_HIV_HXB2 import load_adapter_table



# Globals
VERBOSE = 1
data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/'
data_folder_subsample = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/subsample/'
adapters_table_file = 'adapters_table.dat'
maxreads = 30000

# Quality threshold
phred_min = 20
block_len_min = 100




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

        # Prepare output data structures
        reads1 = []
        reads2 = []

        # Prepare temporary data structures
        read_good = np.zeros(2, bool)
        
        # Directory to read
        dirname = 'adapterID_'+'{:02d}'.format(adaID)+'/'

        # Scroll read files with demultiplexed raw reads
        read1_it = SeqIO.parse(data_folder+dirname+'read1.fastq', 'fastq')
        read2_it = SeqIO.parse(data_folder+dirname+'read2.fastq', 'fastq')
        reads_out = 0
        for seqs12 in izip(read1_it, read2_it):

            # Limit to a small subsample of good reads
            if reads_out >= maxreads:
                break

            # check both reads together
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
                reads1.append(reads_trimmed[0])
                reads2.append(reads_trimmed[1])
                reads_out += 1

                # Keep a log
                if VERBOSE:
                    if not (reads_out % 1000):
                        print reads_out

        # Save output to file
        if reads_out:
            SeqIO.write(reads1, data_folder_subsample+dirname+'read1_filtered_trimmed.fastq', 'fastq')
            SeqIO.write(reads2, data_folder_subsample+dirname+'read2_filtered_trimmed.fastq', 'fastq')
