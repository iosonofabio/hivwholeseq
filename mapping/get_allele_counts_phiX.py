# vim: fdm=indent
'''
author:     Richard Neher, Fabio Zanini
date:       07/10/13
content:    Get the allele counts for phiX (control)
'''
# TODO: test it
import os
import cPickle as pickle
import numpy as np
from Bio import SeqIO
import pysam
import argparse
import get_allele_counts

from mapping.datasets import MiSeq_runs
from mapping.filenames import get_phix_filename, get_mapped_phix_filename, \
        get_allele_counts_phix_filename, get_insert_counts_phix_filename
from mapping.mapping_utils import convert_sam_to_bam



# Globals
alpha = np.array(list('ACGT-N'), 'S1')
read_types = ['read1 f', 'read1 r', 'read2 f', 'read2 r']
maxreads = 5e10
match_len_min = 30
trim_bad_cigars = 3



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Divide HIV reads into fragments')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    miseq_run = args.run
    VERBOSE = args.verbose

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # Read reference
    ref_file = get_phix_filename()
    refseq = SeqIO.read(ref_file, 'fasta')
    ref = np.array(refseq)

    # Open BAM
    bamfilename = get_mapped_phix_filename(data_folder, filtered=True)
    # Try to convert to BAM if needed
    if not os.path.isfile(bamfilename):
        convert_sam_to_bam(bamfilename)

    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        # Get counts
        counts, inserts = get_allele_counts.get_allele_counts(ref, bamfile)
    
    # Save counts and inserts to file
    allele_count_filename = get_allele_counts_phix_filename(data_folder)
    insert_count_filename = get_insert_counts_phix_filename(data_folder)
    np.save(allele_count_filename, counts)
    with open(insert_count_filename, 'w') as f:
        pickle.dump(inserts, f, protocol=-1)
