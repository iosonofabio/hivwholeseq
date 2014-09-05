# vim: fdm=marker
'''
author:     Fabio Zanini
date:       20/03/14
content:    Test script to measure decorrelation (recombination): pick pairs of
            alleles that are strongly correlated at the first time point, and see
            how fast the correlation disappears - if at all.
'''
# Modules
import sys
import os
import argparse
import numpy as np
import pysam
from Bio import SeqIO

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.patients import get_patient
from hivwholeseq.sequencing.filter_mapped_reads import plot_distance_histogram, \
        plot_distance_histogram_sliding_window, get_distance_from_consensus, \
        check_overhanging_reads, trim_bad_cigar
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_mapped_to_initial_filename, get_filter_mapped_init_summary_filename, \
        get_allele_cocounts_filename, get_allele_frequency_trajectories_filename
from hivwholeseq.mapping_utils import convert_sam_to_bam, pair_generator
from hivwholeseq.two_site_statistics import get_coallele_counts_from_file




# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Filter mapped reads')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose

    patient = get_patient(pname)
    samplenames = patient.samples[:2] #FIXME # They are time sorted already
    times = patient.times()

    # If the script is called with no fragment, iterate over all
    # NOTE: This assumes that the fragment has been sequenced. If not, we should
    # take a later time point or similia, but let's not worry about it yet.
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    for fragment in fragments:
        #refseq = SeqIO.read(get_initial_reference_filename(pname, fragment), 'fasta')

        #if VERBOSE >= 1:
        #    print 'Initializing matrix of coallele frequencies'
        #coaft = np.empty((len(samplenames), len(alpha), len(alpha), len(refseq), len(refseq)), float)
        #if VERBOSE >= 1:
        #    print 'Collecting cocounts and normalizing'
        #for i, samplename in enumerate(samplenames):
        #    if VERBOSE >= 2:
        #        print pname, fragment, samplename
        #    cocounts = np.load(get_allele_cocounts_filename(pname, samplename, fragment))
        #    # NOTE: uncovered will get 0 correlation
        #    coaft[i] = 1.0 * cocounts / (cocounts.sum(axis=0).sum(axis=0) + 0.1)
        #    del cocounts

        #if VERBOSE >= 1:
        #    print 'Getting allele frequencies'
        #aft = np.load(get_allele_frequency_trajectories_filename(pname, fragment))[:coaft.shape[0]]
        #if VERBOSE >= 1:
        #    print 'Broadcasting allele frequency product (outer tensor product)'
        #aftbroad = np.einsum('mij...,m...kl->mikjl', aft, np.ones_like(aft))
        #aftbroadT = aftbroad.swapaxes(1, 2).swapaxes(3, 4)
        #if VERBOSE >= 1:
        #    print 'Subtract product of allele frequencies'
        #LD = coaft - (aftbroad * aftbroadT)

        ## Mask LD when allele freqs are too close to 0 or 1
        #LD = np.ma.array(LD)
        #LD[(aftbroad < 1e-2) | (aftbroad > 1 - 1e-2) | \
        #   (aftbroadT < 1e-2) | (aftbroadT > 1 - 1e-2)] = np.ma.masked

        ## Take the correlation coefficient, i.e. divide by sqrt(pi qi pj qj)
        #corr = LD / np.sqrt(aftbroad * (1 - aftbroad) * aftbroadT * (1 - aftbroadT))

        if VERBOSE >= 1:
            print 'Finding initial highly correlated pairs'

        # FIXME: only use this criterion because we have t0 twice
        ind_highLD0 = (corr[0] > 1).nonzero()
        #ind_highLD0 = ((LD[0] > 0.1) & (LD[1] > 0.1) & (aftbroad[0] > 1e-3) & (aftbroad[0] < 0.5)).nonzero()
        pairs_high_LD0 = zip(*ind_highLD0)
        LD_pairs = np.array([LD[:, pair[0], pair[1], pair[2], pair[3]] for pair in pairs_high_LD0]).T
        corr_pairs = np.array([corr[:, pair[0], pair[1], pair[2], pair[3]] for pair in pairs_high_LD0]).T
        print corr_pairs.T
