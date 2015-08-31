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
from itertools import izip
from operator import itemgetter
import numpy as np
from Bio import SeqIO

from hivwholeseq.utils.miseq import alpha
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.one_site_statistics import plot_coverage_trajectories_heatmap, \
        plot_coverage_trajectories_3d, plot_coverage_trajectories



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get coverage trajectories',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+', default=['all'],
                        help='Patients to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the allele frequency trajectories')
    parser.add_argument('--PCR1', type=int, default=1,
                        help='Take only PCR1 samples [0=off, 1=where both available, 2=always]')

    args = parser.parse_args()
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    plot = args.plot
    use_PCR1 = args.PCR1

    patients = load_patients()
    if pnames != ['all']:
        patients = patients.iloc[patients.index.isin(pnames)]

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()

        for fragment in fragments:
            if VERBOSE >= 1:
                print pname, fragment

            covt, ind = patient.get_coverage_trajectories(fragment,
                                                          use_PCR1=use_PCR1)
            samples = patient.samples.iloc[ind]
            times = patient.times[ind]
            ntemplates = samples['n templates']

            if plot is not None:
                import matplotlib.pyplot as plt
                
                if plot in ('2D', '2d', ''):
                    labels = [str(t)+': '+s for (t, s) in izip(times, samples.index)]
                    legtitle = 'Time from transmission [days]: sample name'
                    plot_coverage_trajectories(times, covt,
                                               title='Patient '+pname+', '+fragment,
                                               labels=labels,
                                               legendtitle=legtitle,
                                               VERBOSE=VERBOSE)
                    #plt.tight_layout()
                    #plt.savefig('/ebio/ag-neher/home/fzanini/phd/sequencing/figures/'+\
                    #            'coverage_'+pname+'_'+fragment+'.png')

                elif plot in ('heatmap', 'heat'):
                    plot_coverage_trajectories_heatmap(times, covt,
                                                       title='Patient '+pname+', '+fragment,
                                                       VERBOSE=VERBOSE)

                elif plot in ('3D', '3d'):
                    plot_coverage_trajectories_3d(times, covt,
                                                  title='Patient '+pname+', '+fragment,
                                                  VERBOSE=VERBOSE)

    if plot:
        import matplotlib.pyplot as plt
        plt.ion()
        plt.show()
