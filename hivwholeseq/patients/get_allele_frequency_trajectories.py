#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/11/13
content:    Collect the allele frequencies of all samples to the initial
            consensus by collecting from matrices of single patient samples.
'''
# Modules
import os
import argparse
from operator import itemgetter
import numpy as np
from Bio import SeqIO

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.patients import load_patient
from hivwholeseq.patients.filenames import get_initial_consensus_filename, \
        get_mapped_to_initial_filename, get_allele_frequency_trajectories_filename, \
        get_allele_count_trajectories_filename
from hivwholeseq.patients.one_site_statistics import \
        plot_allele_frequency_trajectories_from_counts as plot_nus_from_act
from hivwholeseq.patients.one_site_statistics import \
        plot_allele_frequency_trajectories_from_counts_3d as plot_nus_from_act_3d
from hivwholeseq.patients.one_site_statistics import get_allele_count_trajectories
from hivwholeseq.fork_cluster import fork_get_allele_frequency_trajectory as fork_self



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get allele frequency trajectories')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    #parser.add_argument('--submit', action='store_true',
    #                    help='Execute the script in parallel on the cluster')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the allele frequency trajectories')
    parser.add_argument('--logit', action='store_true',
                        help='use logit scale (log(x/(1-x)) in the plots')
    parser.add_argument('--PCR1', type=int, default=1,
                        help='Take only PCR1 samples [0=off, 1=where both available, 2=always]')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose
    plot = args.plot
    #submit = args.submit
    use_PCR1 = args.PCR1
    use_logit = args.logit

    patient = load_patient(pname)
    patient.discard_nonsequenced_samples()
    samplenames = patient.samples.index

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    # Iterate over samples and fragments
    for fragment in fragments:

        ## Submit to the cluster self if requested
        #if submit:
        #    fork_self(pname, fragment, VERBOSE=VERBOSE)
        #    continue

        if VERBOSE >= 1:
            print fragment

        act_filename = get_allele_count_trajectories_filename(pname, fragment)
        aft_filename = get_allele_frequency_trajectories_filename(pname, fragment)

        # Collect allele counts from patient samples, and return only positive hits
        # sns contains sample names and PCR types
        (sns, act) = get_allele_count_trajectories(pname, samplenames, fragment,
                                                   use_PCR1=use_PCR1, VERBOSE=VERBOSE)

        ind = [i for i, (_, sample) in enumerate(patient.samples.iterrows())
               if sample.name in map(itemgetter(0), sns)]
        samples = patient.samples.iloc[ind]
        times = (samples.date - patient.transmission_date) / np.timedelta64(1, 'D')

        # FIXME: use masked arrays?
        aft = 1.0 * act / act.sum(axis=0)
        aft[np.isnan(aft)] = 0
        aft[(aft < 1e-5) | (aft > 1)] = 0

        if plot is not None:
            import matplotlib.pyplot as plt

            if plot in ('2D', '2d', ''):

                # FIXME: the number of molecules to PCR depends on the number of
                # fragments for that particular experiment... integrate Lina's table!
                # Note: this refers to the TOTAL # of templates, i.e. the factor 2x for
                # the two parallel RT-PCR reactions
                ntemplates = samples['viral load'] * 0.4 / 12 * 2

                plot_nus_from_act(times, act,
                                  title='Patient '+pname+', '+fragment,
                                  VERBOSE=VERBOSE, logit=use_logit,
                                  ntemplates=ntemplates,
                                  threshold=0.9)

            elif plot in ('3D', '3d'):
                plot_nus_from_act_3d(times, act,
                                     title='Patient '+pname+', '+fragment,
                                     VERBOSE=VERBOSE, logit=use_logit,
                                     threshold=0.9)

    if plot is not None:
        plt.tight_layout()
        plt.ion()
        plt.show()
