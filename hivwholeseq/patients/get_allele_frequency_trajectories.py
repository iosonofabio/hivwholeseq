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
from hivwholeseq.patients.patients import load_patient, convert_date_deltas_to_float
from hivwholeseq.patients.filenames import get_allele_frequency_trajectories_filename, \
        get_allele_count_trajectories_filename
from hivwholeseq.patients.one_site_statistics import \
        plot_allele_frequency_trajectories as plot_aft
from hivwholeseq.patients.one_site_statistics import \
        plot_allele_frequency_trajectories_3d as plot_aft_3d
from hivwholeseq.patients.one_site_statistics import get_allele_count_trajectories
from hivwholeseq.fork_cluster import fork_get_allele_frequency_trajectory as fork_self



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get allele frequency trajectories',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 genomewide)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the allele frequency trajectories')
    parser.add_argument('--logit', action='store_true',
                        help='use logit scale (log(x/(1-x)) in the plots')
    parser.add_argument('--PCR1', type=int, default=1,
                        help='Take only PCR1 samples [0=off, 1=where both available, 2=always]')
    parser.add_argument('--threshold', type=float, default=0.9,
                        help='Minimal frequency plotted')
    parser.add_argument('--cov-min', type=float, default=200,
                        help='Minimal coverage (lower sites are masked)')
    parser.add_argument('--depth-min', type=float, default=200,
                        help='Minimal depth (lower time points are hidden)')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose
    plot = args.plot
    use_PCR1 = args.PCR1
    use_logit = args.logit
    threshold = args.threshold
    cov_min = args.cov_min
    depth_min = args.depth_min

    patient = load_patient(pname)
    patient.discard_nonsequenced_samples()
    samplenames = patient.samples.index

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    for fragment in fragments:
        if VERBOSE >= 1:
            print fragment

        # Collect allele counts from patient samples, and return only positive hits
        # sns contains sample names and PCR types
        (aft, ind) = patient.get_allele_frequency_trajectories(fragment,
                                                           use_PCR1=use_PCR1,
                                                           cov_min=cov_min,
                                                           depth_min=depth_min,
                                                           VERBOSE=VERBOSE)
        times = patient.times[ind]
        ntemplates = patient.n_templates[ind]

        if plot is not None:
            import matplotlib.pyplot as plt

            if plot in ('2D', '2d', ''):
                (fig, ax) = plot_aft(times, aft,
                                     title='Patient '+pname+', '+fragment,
                                     VERBOSE=VERBOSE,
                                     ntemplates=ntemplates,
                                     logit=use_logit,
                                     threshold=threshold)

            elif plot in ('3D', '3d'):
                (fig, ax) = plot_aft_3d(times, aft,
                                        title='Patient '+pname+', '+fragment,
                                        VERBOSE=VERBOSE,
                                        logit=use_logit,
                                        ntemplates=ntemplates,
                                        threshold=threshold)

    if plot is not None:
        plt.tight_layout()
        plt.ion()
        plt.show()
