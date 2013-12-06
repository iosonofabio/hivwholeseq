# vim: fdm=marker
'''
author:     Fabio Zanini
date:       25/08/13
content:    Test script for extracting linkage information from the reads.
'''
# Modules
import os
import sys
import argparse
import cPickle as pickle
from collections import defaultdict
from itertools import izip
import pysam
import numpy as np
from Bio import SeqIO

from hivwholeseq.filenames import get_last_reference, get_last_mapped, get_mutations_file



# Globals
VERBOSE = 1

# FIXME
from hivwholeseq.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']




# Functions
def load_mutations(data_folder, adaID):
    '''Load the pickle file with the mutations'''
    mut_file = get_mutations_file(data_folder, adaID)
    from time import time
    t0 = time()
    with open(mut_file, 'r') as f:
        muts_all = pickle.load(f)
    t1 = time()
    print int(t1 - t0)
    return muts_all


def filter_only_mutations_seen_x_times(muts_all, times=4):
    '''Extract mutations that have been seen at least x times'''
    from collections import Counter

    # 1. Merge mutations from all reads
    muts_flat = []
    for mut in muts_all:
        # mut[1] is the fragment, m[0] the site, m[1] the derived allele
        muts_flat.extend(((mut[1], m[0], m[1]) for m in mut[2]))

    # 2. Count and filter them
    allele_counts = Counter(muts_flat)
    alleles_good = {m: count for m, count in allele_counts.iteritems() if count >= times}
    return alleles_good



# Script
if __name__ == '__main__':

#    # Input arguments
#    parser = argparse.ArgumentParser(description='Extract linkage information')
#    parser.add_argument('--adaID', metavar='00', type=int, required=True,
#                        help='Adapter ID sample to analyze')
#    args = parser.parse_args()
#    adaID = args.adaID
#
#    # Load the mutations
#    muts_all = load_mutations(data_folder, adaID)

    # Pick the mutations seen at least 4 times (the rest is error)
    from time import time
    t0 = time()
    mutations_good = filter_only_mutations_seen_x_times(muts_all, 4)
    t1 = time()
    print int(t1 - t0)

    # Sort them by frequency
    muts, counts = zip(*mutations_good.iteritems())
    ind = np.argsort(counts)[::-1]
    muts = np.array(muts, dtype=list)[ind]
    counts = np.array(counts)[ind]


    # TODO: exclude known high-error sites (Richard should have written the
    # error filters somewhere?)
