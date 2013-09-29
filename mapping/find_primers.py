# vim: fdm=indent
'''
author:     Fabio Zanini
date:       25/09/13
content:    Reconstruct the primers coordinates (approx.) from the edges of the
            fragments.
'''
# Modules
import os
import sys
import argparse
import cPickle as pickle
from operator import itemgetter
from collections import defaultdict
from itertools import izip
import pysam
import numpy as np
from Bio import SeqIO


# Horizontal import of modules from this folder
from mapping.adapter_info import load_adapter_table, foldername_adapter
from mapping.miseq import alpha, read_types
from mapping.filenames import get_consensus_filename, get_mapped_filename
from mapping.mapping_utils import get_ind_good_cigars, convert_sam_to_bam,\
        pair_generator, get_range_good_cigars



# Globals
# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Extract linkage information')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--verbose', type=int, default=0,
                        help=('Verbosity level [0-3]'))
    args = parser.parse_args()
    adaIDs = args.adaIDs
    VERBOSE = args.verbose

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = load_adapter_table(data_folder)['ID']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # Select fragment and primers
    fragment = 'F3'
    # Look for the F3 rev primer (already reversed)
    primer_old = 'GATTGTGTGGCAAGTAGACAGG'
    primer_new = 'TATGGAAAACAGATGGCAGGTG'

    # Iterate over all requested samples
    for adaID in adaIDs:

        # Read reference (fragmented)
        reffilename = get_consensus_filename(data_folder, adaID, fragment)
        refseq = SeqIO.read(reffilename, 'fasta')
        ref = np.array(refseq)

        # read file
        bamfilename = get_mapped_filename(data_folder, adaID, fragment,
                                          type='bam', filtered=True)

        if not os.path.isfile(bamfilename):
            convert_sam_to_bam(bamfilename)
        bamfile = pysam.Samfile(bamfilename, 'rb')

        # Get the coverage for reads which have long insert sizes
        # (to be sure about their identity)
        cov_new = 0
        cov_old = 0
        for i_pairs, reads in enumerate(pair_generator(bamfile)):
            if i_pairs > 5000000:
                break

            if reads[0].isize < 300:
                continue

            for read in reads:
                if read.seq.find(primer_new) != -1:
                    cov_new += 1
                if read.seq.find(primer_old) != -1:
                    cov_old += 1

        print 'old:', cov_old, 'new:', cov_new

        bamfile.close()

        # Get coverage and see
        covfn = get_coverage_filename(data_folder, adaID, fragment)
        cov = np.load(covfn)

        import matplotlib.pyplot as plt
        import matplotlib.cm as cm

        for js, read_type in enumerate(read_types):
            plt.plot(np.arange(cov.shape[1]), cov[js], lw=2,
                     c=cm.jet(int(255.0 * js / len(read_types))))

        plt.xlabel('Position [bases]')
        plt.title(str(adaID)+' '+fragment)
        plt.ylabel('Coverage')

        plt.ion()
        plt.show()

