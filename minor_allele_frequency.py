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



# Globals
VERBOSE = 1
alpha = np.array(list('ACGT-N'), 'S1')
read_types = ['read1 f', 'read1 r', 'read2 f', 'read2 r']
ref_filename = 'consensus_filtered_trimmed.fasta'
allele_count_filename = 'allele_counts.npy'



# Script
if __name__ == '__main__':

    # Get data
    adaID = 02
    data_folder = ('/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/'+
                   'subsample/'+
                   'adapterID_'+'{:02d}'.format(adaID)+'/')
    ref_file = data_folder+ref_filename

    # Read reference
    # Read reference
    if os.path.isfile(data_folder+ref_filename): ref_file = data_folder+ref_filename
    else: ref_file = '/'.join(data_folder.split('/')[:-2]+['subsample/']+data_folder.split('/')[-2:])+ref_filename
    refseq = SeqIO.read(ref_file, 'fasta')[1000:8500]
    ref = np.array(refseq)

    # Get counts
    # The first dimension is read type, the second alphabet, the third position
    counts = np.load(data_folder+allele_count_filename)[:, :, 1000:8500]

    # Get minor allele frequencies
    coverage = counts.sum(axis=1)
    major = counts.argmax(axis=1)
    nu_major = 1.0 * counts.max(axis=1) / (coverage + 0.0000001)
    # 1. Main minor allele
    counts_minor = counts.copy()
    for pos in xrange(counts_minor.shape[-1]):
        for js in xrange(counts_minor.shape[0]):
            counts_minor[js, major[js, pos], pos] = -1
    minor = counts_minor.argmax(axis=1)
    nu_minor = 1.0 * counts_minor.max(axis=1) / (coverage + 0.0000001)
    # 2. Second minor allele
    counts_minor2 = counts_minor.copy()
    for pos in xrange(counts_minor.shape[-1]):
        for js in xrange(counts_minor.shape[0]):
            counts_minor2[js, minor[js, pos], pos] = -1
    minor2 = counts_minor2.argmax(axis=1)
    nu_minor2 = 1.0 * counts_minor2.max(axis=1) / (coverage + 0.0000001)

    # Filter out sequencing errors
    errors_seq = np.zeros(coverage.shape[1], bool)

    # Sequencing errors tend to have this pattern: they are seen in both forward
    # reads, but not reverse reads (or vice versa)
    nu_mav = nu_minor.mean(axis=0)
    corr_sense = (nu_minor[0] - nu_mav) * (nu_minor[2] - nu_mav) + (nu_minor[1] - nu_mav) * (nu_minor[3] - nu_mav)
    corr_antis = (nu_minor[0] - nu_mav) * (nu_minor[1] - nu_mav) + (nu_minor[2] - nu_mav) * (nu_minor[3] - nu_mav)
    errors_seq |= (corr_sense > 0) & (corr_antis < 0) & (nu_mav > 0.01)
    
    errors_seq |= (nu_minor > 0.01).sum(axis=0) < 2 # Ignore if they are not seen in at least 2 read types

    ## Tag as errors if their frequencies in different read types are too different
    #for pos in xrange(len(errors_seq)):
    #    nu_seen = nu_minor[nu_minor[:, pos] > 0.01, pos]
    #    if len(nu_seen) and (nu_seen.std() / nu_seen.mean() > 0.1):
    #        errors_seq[pos] = True

    # Tag as errors if the coverage is real small
    errors_seq |= (coverage < 50).any(axis=0)

    #plt.scatter(corr_sense, corr_antis, s=30)
    #plt.plot([-0.15, 0.15], [0] * 2, lw=1.5, ls='--', c='k')
    #plt.plot([0] * 2, [-0.15, 0.15], lw=1.5, ls='--', c='k')


    # Plot?
    fig, axs = plt.subplots(2, 1, figsize=(14, 10))
    for i, (nu_m, read_type) in enumerate(izip(nu_minor, read_types)):
        y = nu_m.copy()
        y[errors_seq] = 0
        axs[0].plot(y + 0.0000001, lw=1.5)
        h = np.histogram(y, bins=np.logspace(-3, 0, 100), density=False)
        axs[1].plot(np.sqrt(h[1][1:] * h[1][:-1]), h[0], lw=2, label=read_type)

    #axs[0].set_ylim(0.05, 1.05)
    axs[0].set_yscale('log')
    axs[1].set_xscale('log')
    axs[1].set_xlim(0.95e-3, 1.05)
    axs[1].set_ylim(-0.15, axs[1].get_ylim()[1])
    axs[0].set_xlabel('HIV genome [b.p.]')
    axs[0].set_ylabel('Minor allele frequency')
    axs[1].set_xlabel('Minor allele frequency')
    axs[1].set_ylabel('# alleles (should be flat @ zero)')
    axs[1].legend(loc=1)
    plt.ion()
    plt.show()
