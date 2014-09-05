# vim: fdm=marker
'''
author:     Fabio Zanini
date:       28/03/14
content:    Test script for measuring recombination.
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
    samplenames = patient.samples
    times = patient.times()

    # Keep PCR2 only if PCR1 is absent
    criterion = lambda x: ('PCR1' in x[1]) or ((times == times[x[0]]).sum() == 1)
    ind = np.nonzero(map(criterion, enumerate(samplenames)))[0]
    times = times[ind]

    # If the script is called with no fragment, iterate over all
    # NOTE: This assumes that the fragment has been sequenced. If not, we should
    # take a later time point or similia, but let's not worry about it yet.
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    for fragment in fragments:
        refseq = SeqIO.read(get_initial_reference_filename(pname, fragment), 'fasta')

        if VERBOSE >= 1:
            print 'Collecting initial cocounts and cofrequencies'
        cocounts = np.ma.array(np.load(get_allele_cocounts_filename(pname, samplenames[0], fragment)))
        coaf = 1.0 * cocounts / (cocounts.sum(axis=0).sum(axis=0))
        coaf[cocounts < 100] = np.ma.masked

        if VERBOSE >= 1:
            print 'Getting allele frequencies'
        aft = np.load(get_allele_frequency_trajectories_filename(pname, fragment))[ind]
        aftbroad = np.einsum('mij...,m...kl->mikjl', aft, np.ones_like(aft))
        aftbroadT = aftbroad.swapaxes(1, 2).swapaxes(3, 4)
        if VERBOSE >= 1:
            print 'Calculate LD'
        LD = np.ma.array(coaf - (aftbroad[0] * aftbroadT[0]))

        if VERBOSE:
            print 'Getting distances'
        dist = np.tile(np.arange(aft.shape[-1]), (aft.shape[-1], 1))
        dist = np.abs(dist - dist.T)
        distbroad = np.tile(dist, (aft.shape[1], aft.shape[1], 1, 1))

        if VERBOSE:
            print 'Selecting pairs'
        mask = (aftbroad[0] < 1e-2) | (aftbroad[0] > 1 - 1e-2) | \
               (aftbroadT[0] < 1e-2) | (aftbroadT[0] > 1 - 1e-2) | \
               (coaf < 1e-3) | (distbroad == 0)
        LD[mask] = np.ma.masked

        if VERBOSE:
            print 'Calculating LD sign'
        sign = np.ma.masked_all(LD.shape)
        sign[LD > 0.001] = 1
        sign[LD < -0.001] = -1

        if VERBOSE:
            print 'Calculate linkage statistics: sgn(LD0) * Dnu_i * Dnu_j'
        st = sign * (aftbroad[1:] - aftbroad[:-1]) * (aftbroadT[1:] - aftbroadT[:-1])
        
        if VERBOSE:
            print 'Classifying by distance'
        binwidth = 30
        dbins = distbroad // binwidth
        ds = binwidth // 2 + binwidth * np.arange(5)
        stav = [st[:, dbins == i].mean() for i in xrange(5)]
        print ds; print map('{:2.3e}'.format, stav)
