# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/08/13
content:    Check sequencing/PCR erros in NL4-3 (they appear as minor variants).
'''
# Modules
import os
import sys
from collections import defaultdict
from itertools import izip
import pysam
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from get_allele_frequencies import get_allele_counts



# Globals
VERBOSE = 1
alpha = np.array(list('ACGT-N'), 'S1')
read_types = ['read1 f', 'read1 r', 'read2 f', 'read2 r']
maxreads = 20000
ref_filename = 'consensus_filtered_trimmed.fasta'
bam_filename = 'mapped_to_self_filtered_trimmed.bam'



# Script
if __name__ == '__main__':

    # Get data
    data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/subsample/adapterID_02/'
    ref_file = data_folder+ref_filename
    bamfilename = data_folder+bam_filename

    # Read reference
    refseq = SeqIO.read(ref_file, 'fasta')
    ref = np.array(refseq)

    # Open BAM
    bamfile = pysam.Samfile(bamfilename, 'rb')
    
    # Get counts
    # The first dimension is read type, the second alphabet, the third position
    counts, inserts = get_allele_counts(ref, bamfile)

    # Get minor allele frequencies
    counts_minor = counts.copy()
    coverage = counts.sum(axis=1)
    nu_major = counts.max(axis=1) / (coverage + 0.000001)
