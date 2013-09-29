# vim: fdm=marker
'''
author:     Fabio Zanini
date:       16/09/13
content:    Study the read length distribution at the end of the mapping pipeline.
'''
# Modules
import os
import sys
import cPickle as pickle
import argparse
import re
from collections import defaultdict
from itertools import izip
import pysam
import numpy as np
from Bio import SeqIO
import matplotlib.cm as cm

# Matplotlib parameters
import matplotlib
params = {'axes.labelsize': 20, 
          'text.fontsize': 20,
          'legend.fontsize': 18,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16,
          'text.usetex': False}
matplotlib.rcParams.update(params)
import matplotlib.pyplot as plt


# Globals
# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']

from mapping.adapter_info import load_adapter_table
from mapping.miseq import alpha, read_types
from mapping.filenames import get_mapped_filename, get_allele_counts_filename, \
        get_insert_counts_filename, get_coverage_filename, get_consensus_filename
from mapping.mapping_utils import get_fragment_list, convert_sam_to_bam



def get_read_lengths(data_folder, adaID, fragment, subsample=False, VERBOSE=0):
    '''Get the read lengths'''

    # Lengths from 1 to 250
    lengths = np.zeros((len(read_types), 250), int)

    # Open BAM file
    # Note: the reads should already be filtered of unmapped stuff at this point
    bamfilename = get_mapped_filename(data_folder, adaID, fragment, type='bam',
                                      filtered=True)
    if not os.path.isfile(bamfilename):
        convert_sam_to_bam(bamfilename)
    with pysam.Samfile(bamfilename, 'rb') as bamfile:

        # Iterate over single reads (no linkage info needed)
        for i, read in enumerate(bamfile):
        
            # Print output
            if (VERBOSE >= 3) and (not ((i +1) % 10000)):
                print (i+1)
        
            # Divide by read 1/2 and forward/reverse
            js = 2 * read.is_read2 + read.is_reverse

            # Increment counter
            lengths[js, read.rlen - 1] += 1

            # Note: we do not delve into CIGARs because the reads are trimmed

    return lengths



# Script
if __name__ == '__main__':


    # Input arguments
    parser = argparse.ArgumentParser(description='Study minor allele frequency')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--subsample', action='store_true',
                        help='Apply only to a subsample of the reads')

    args = parser.parse_args()
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    subsample = args.subsample

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = load_adapter_table(data_folder)['ID']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Iterate over all requested samples
    lengths_all = {}
    for adaID in adaIDs:
        for fragment in fragments:
            lengths_all[(adaID, fragment)] = length = \
                    get_read_lengths(data_folder, adaID, fragment,
                                     subsample=subsample, VERBOSE=VERBOSE)

            # Plot it
            fig, ax = plt.subplots(1, 1)
            for irt, read_type in enumerate(read_types):
                color = cm.jet(int(255.0 * irt / len(read_types)))
                ax.plot(np.arange(1, 251), length[irt], lw=1.5, c=color)
                ax.scatter(np.arange(1, 251), length[irt], s=50, c=color)
            ax.set_xlabel('Read length [bases]')
            ax.set_ylabel('#')
            ax.set_title('{:02d}'.format(adaID)+', '+fragment)
            ax.set_yscale('log')
            
            plt.tight_layout()
            plt.show()


