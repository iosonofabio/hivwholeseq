# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/08/13
content:    Extract the frequency of minor alleles.
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
from map_HIV_HXB2 import load_adapter_table


# Globals
VERBOSE = 1
alpha = np.array(list('ACGT-N'), 'S1')
read_types = ['read1 f', 'read1 r', 'read2 f', 'read2 r']
ref_filename = 'consensus_filtered_trimmed.fasta'
allele_count_filename = 'allele_counts.npy'



# Functions
def get_minor_allele_counts(counts, n_minor=0):
    '''Extract the minor allele counts
    
    Parameters:
       - counts: the complete allele counts
       - n_minor: the top minor allele or the second-top?
    '''
    
    # 0. Major allele
    major = counts.argmax(axis=1)
    # 1. Main minor allele
    counts_minor = counts.copy()
    for pos in xrange(counts_minor.shape[-1]):
        for js in xrange(counts_minor.shape[0]):
            counts_minor[js, major[js, pos], pos] = -1
    minor = counts_minor.argmax(axis=1)
    counts_minor = counts_minor.max(axis=1)
    if not n_minor:
        return minor, counts_minor
    elif n_minor == 1:
        # 2. Second minor allele
        for pos in xrange(counts_minor.shape[-1]):
            for js in xrange(counts_minor.shape[0]):
                counts_minor[js, minor[js, pos], pos] = -1
        minor2 = counts_minor.argmax(axis=1)
        counts_minor2 = counts_minor.max(axis=1)
        return minor2, counts_minor2
    else:
        raise NotImplemented('4+th allele not implemented')




# Script
if __name__ == '__main__':

    # Input args
    args = sys.argv
    if len(args) < 2:
        raise ValueError('Please select an adapterID')
    adaID = int(args[1])

    # Get data
    data_folder = ('/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/')
    data_adaID_folder = data_folder+'adapterID_'+'{:02d}'.format(adaID)+'/'

    # Get adapter table
    adapter_table = load_adapter_table(data_folder)
    sample_name = adapter_table['sample'][adapter_table['ID'] == adaID][0]

    ## Read reference
    #ref_file = data_adaID_folder+ref_filename
    #if os.path.isfile(data_adaID_folder+ref_filename): ref_file = data_adaID_folder+ref_filename
    #else: ref_file = '/'.join(data_adaID_folder.split('/')[:-2]+['subsample/']+data_adaID_folder.split('/')[-2:])+ref_filename
    #refseq = SeqIO.read(ref_file, 'fasta')
    #ref = np.array(refseq)

    # Get counts
    # The first dimension is read type, the second alphabet, the third position
    counts = np.load(data_adaID_folder+allele_count_filename)

    # Define sequenced region
    sequenced_region = [(counts.sum(axis=1)[2] > 9000).nonzero()[0][0] - 1,
                        (counts.sum(axis=1)[3] > 9000).nonzero()[0][-1] + 2]
    counts = counts[:, :, sequenced_region[0]: sequenced_region[1]]

    # Get minor allele frequencies
    coverage = counts.sum(axis=1)
    major = counts.argmax(axis=1)
    nu_major = 1.0 * counts.max(axis=1) / (coverage + 0.0000001)
    # Top minor allele
    minor, counts_minor = get_minor_allele_counts(counts)
    nu_minor = 1.0 * counts_minor / (coverage + 0.0000001)

    # Filter out sequencing errors
    errors_seq = np.zeros(coverage.shape[1], bool)

    ## Exclude SNP calls in regions of low coverage
    #errors_seq |= (coverage < 1000).all(axis=0)

    ## Sequencing errors are characterized by not being present in all four read
    ## types. We check all the types with enough coverage and test for homogeneously
    ## distributed counts
    #for pos in (-errors_seq).nonzero()[0]:
    #    read_types_cov = coverage[:, pos] > 1000
    #    count_pos = counts_minor[read_types_cov, pos]
    #    if count_pos.any():
    #        cov = coverage[read_types_cov, pos]
            

    ## Sequencing errors tend to have this pattern: they are seen in both forward
    ## reads, but not reverse reads (or vice versa)
    #nu_mav = nu_minor.mean(axis=0)
    #corr_sense = (nu_minor[0] - nu_mav) * (nu_minor[2] - nu_mav) + (nu_minor[1] - nu_mav) * (nu_minor[3] - nu_mav)
    #corr_antis = (nu_minor[0] - nu_mav) * (nu_minor[1] - nu_mav) + (nu_minor[2] - nu_mav) * (nu_minor[3] - nu_mav)
    #errors_seq |= (corr_sense > 0) & (corr_antis < 0) & (nu_mav > 0.001)
    #
    #errors_seq |= (nu_minor > 0.01).sum(axis=0) < 2 # Ignore if they are not seen in at least 2 read types

    ### Tag as errors if their frequencies in different read types are too different
    ##for pos in xrange(len(errors_seq)):
    ##    nu_seen = nu_minor[nu_minor[:, pos] > 0.01, pos]
    ##    if len(nu_seen) and (nu_seen.std() / nu_seen.mean() > 0.1):
    ##        errors_seq[pos] = True

    ## Tag as errors if the coverage is real small
    #errors_seq |= (coverage < 3000).any(axis=0)

    #plt.scatter(corr_sense, corr_antis, s=30)
    #plt.plot([-0.15, 0.15], [0] * 2, lw=1.5, ls='--', c='k')
    #plt.plot([0] * 2, [-0.15, 0.15], lw=1.5, ls='--', c='k')

    # Plot?
    # Set special stuff for special plots
    if adaID == 16: density = True
    else: density = False

    # Plot data
    import matplotlib.cm as cm
    fig, axs = plt.subplots(3, 1, figsize=(14, 14))
    for i in xrange(len(read_types)):
        axs[0].plot(sequenced_region[0] + np.arange(coverage.shape[1]), coverage[i],
                    c=cm.jet(int(255.0 * i / len(read_types))))
        read_type = read_types[i]
        y = nu_minor[i].copy()
        #y[errors_seq] = 0
        axs[1].plot(sequenced_region[0] + np.arange(coverage.shape[1]), y + 0.0000001,
                    lw=1.5, c=cm.jet(int(255.0 * i / len(read_types))))
        h = np.histogram(y, bins=np.logspace(-3, 0, 100), density=density)
        axs[2].plot(np.sqrt(h[1][1:] * h[1][:-1]), h[0], lw=2, label=read_type, c=cm.jet(int(255.0 * i / len(read_types))))

    axs[0].set_ylabel('Coverage [# reads]')
    axs[1].set_ylim(1e-5, 1e0)
    axs[1].set_yscale('log')
    axs[2].set_xscale('log')
    axs[2].set_xlim(0.95e-3, 1.05)
    axs[2].set_ylim(-0.15, axs[2].get_ylim()[1])
    axs[1].set_xlabel('HIV genome [b.p.]')
    axs[1].set_ylabel('Minor allele frequency')
    axs[2].set_xlabel('Minor allele frequency')
    axs[2].set_ylabel('# alleles')
    axs[2].legend(loc=1)
    axs[0].set_title('Sample: '+sample_name)

    # Set MORE special stuff
    if adaID == 16:
        x = np.logspace(-3, 0, 1000)
        axs[2].plot(x, 1e-3 / x**2, c='k', lw=2, ls='--')
        axs[2].annotate(r'$p(\nu) \sim 1 / \nu^2$', (9e-2, 0.09), (8e-3, 0.04),
                        arrowprops=dict(arrowstyle="->"),
                        fontsize=20)
        axs[2].set_ylim(1e-2, 1e3)
        axs[2].set_yscale('log')
    if adaID == 18:
        axs[1].plot([0, 10000], [2e-1] * 2, lw=1.5, ls='--', c='k', alpha=0.3)
        axs[2].plot([2e-1] * 2, [0, 380], lw=1.5, ls='--', c='k', alpha=0.3)
    if adaID == 19:
        axs[1].plot([0, 10000], [6e-2] * 2, lw=1.5, ls='--', c='k', alpha=0.3)
        axs[1].plot([0, 10000], [1e-2] * 2, lw=1.5, ls='--', c='k', alpha=0.3)
        axs[2].plot([6e-2] * 2, [0, 380], lw=1.5, ls='--', c='k', alpha=0.3)
        axs[2].plot([1e-2] * 2, [0, 380], lw=1.5, ls='--', c='k', alpha=0.3)

    plt.tight_layout()

    plt.ion()
    plt.show()
