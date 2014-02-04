# vim: fdm=marker
'''
author:     Fabio Zanini
date:       03/02/14
content:    Classify premapped reads according to various criteria.
'''
# Modules
import os
import argparse
from itertools import izip
from Bio import SeqIO
import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from hivwholeseq_pacbio.samples import samples
from hivwholeseq_pacbio.datasets import data_folder_dict
from hivwholeseq_pacbio.filenames import get_premapped_file, \
        get_reference_premap_filename


def classify_reads(bamfilename, reflen, maxreads=-1, VERBOSE=0):
    '''Classify reads using various mapping criteria'''

    from collections import defaultdict
    cla = defaultdict(int)

    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        for i, read in enumerate(bamfile):

            if VERBOSE >= 2:
                if not ((i+1) % 1000):
                    print (i+1)
        
            if i == maxreads:
                break
        
            if read.is_unmapped:
                if VERBOSE >= 3:
                    print 'Unmapped'
                cla['unmapped'] += 1
                continue
            
            cla['mapped'] += 1
            
    return cla



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get coverage of PacBio reads')
    parser.add_argument('--run', default='Upp23',
                        help='PacBio run to analyze (e.g. Upp23)')
    parser.add_argument('--sample', required=True,
                        help='Sample to analyze (e.g. S1)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to map')

    args = parser.parse_args()
    seq_run = args.run
    samplename = args.sample
    VERBOSE = args.verbose
    maxreads = args.maxreads

    # Specify the dataset
    data_folder = data_folder_dict[seq_run]
    sample = samples.set_index('name').loc[samplename]

    # Get NL4-3 reference
    refseq = SeqIO.read(get_reference_premap_filename(data_folder, samplename), 'fasta')

    # Get coverage
    bamfilename = get_premapped_file(data_folder, samplename)
    classification = classify_reads(bamfilename, len(refseq), maxreads=maxreads,
                                    VERBOSE=VERBOSE)

    print seq_run, samplename
    print '\t'.join(['Category', '# reads'])
    print '-' * 50
    for (key, value) in classification.iteritems():
        print '\t'.join(['{:15s}'.format(key), '{:1.1e}'.format(value)])
