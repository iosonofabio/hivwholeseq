# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/08/13
content:    Check the coverage of mapped reads on the HIV genome.
'''
# Modules
import os
import sys
import re
import cPickle as pickle
import argparse
from collections import defaultdict, Counter
from itertools import izip
from operator import itemgetter
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt



# Globals
VERBOSE = 1

# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']

from mapping.adapter_info import load_adapter_table
from mapping.miseq import alpha, read_types
from mapping.filenames import get_last_reference, get_coverage_filename
from mapping.mapping_utils import get_fragment_list





# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Extract linkage information')
    parser.add_argument('--adaID', metavar='00', type=int, required=True,
                        help='Adapter ID sample to analyze')
    args = parser.parse_args()
    adaID = args.adaID

    # Fragments
    chromosomes = get_fragment_list(data_folder, adaID)

    # Read reference (fragmented)
    refseqs_raw = list(SeqIO.parse(get_last_reference(data_folder, adaID, ext=True),
                                   'fasta'))
    # Sort according to the chromosomal ordering
    refseqs = []
    for chromosome in chromosomes:
        for seq in refseqs_raw:
            if chromosome == seq.id:
                refseqs.append(seq)
                break
    refs = [np.array(refseq) for refseq in refseqs]

    # Load allele counts
    with open(get_coverage_filename(data_folder, adaID), 'r') as f:
        coverage = pickle.load(f)

    # Argsort chromosomes by name
    ind_chrom = sorted(range(len(chromosomes)), key = chromosomes.__getitem__)

    # Plot the coverage for each fragment
    fig, axs = plt.subplots(2, 3, figsize=(21, 12))
    axs = axs.ravel()
    for ii, i in enumerate(ind_chrom):
        cov = coverage[i]
        ref = refs[i]

        plt.sca(axs[ii])
        for irt, (read_type, co) in enumerate(izip(read_types, cov)):
            plt.plot(co, lw=1.5, c=cm.jet(int(255.0 * irt / len(read_types))),
                     label=read_type)
        if ii > 2:
            plt.xlabel('Position in fragment [b.p.]')
        if ii in [0, 3]:
            plt.ylabel('Coverage')
        plt.title(re.sub('_', ' ', chromosomes[i]))

    fig.suptitle('adapterID '+'{:02d}'.format(adaID), fontsize=20)

    plt.tight_layout(rect=(0, 0, 1, 0.95))

    plt.ion()
    plt.show()
