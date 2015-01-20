# vim: fdm=marker
'''
author:     Fabio Zanini, Richard Neher
date:       19/01/15
content:    Estimate the number of template molecules to PCR, fragment by fragment.
'''
# Modules
import os
import sys
import argparse
from itertools import izip
from collections import defaultdict
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import hivwholeseq.plot_utils

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.filenames import get_ntemplates_by_fragment_filename



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Estimate number of templates from consecutive samples',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='+',
                        default=['F'+str(i) for i in xrange(1, 7)],
                        help='Fragments to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save alignment to file')
    parser.add_argument('--plot', action='store_true',
                        help='Plot changes')

    args = parser.parse_args()
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    use_save = args.save
    use_plot = args.plot

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    data = defaultdict(dict)
    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        for ifr, fragment in enumerate(fragments):
            if VERBOSE >= 1:
                print pname, fragment

            aft, ind = patient.get_allele_frequency_trajectories(fragment,
                                                                 cov_min=100,
                                                                 VERBOSE=VERBOSE)

            times = patient.times[ind]

            for i in xrange(len(times) -1):
                samplename1 = patient.samples.iloc[ind[i]].name
                samplename2 = patient.samples.iloc[ind[i+1]].name

                af1 = aft[i]
                af2 = aft[i+1]
                dt = times[i+1] - times[i]

                if VERBOSE >= 2:
                    print samplename1, 'to', samplename2, '('+str(int(dt))+' days)'

                afmin = 1e-2
                indaf = (af1 >= afmin) & (af1 <= 1 - afmin) & (af2 >= afmin) & (af2 <= 1 - afmin)

                mea = 0.5 * (af1[indaf] + af2[indaf])
                var = (0.5 * (af1[indaf] - af2[indaf]))**2

                eta = (var / (mea * (1 - mea))).mean()

                key = (pname, fragment, ind[i], ind[i+1])
                data['eta'][key] = eta
                data['dt'][key] = dt

                continue

                if use_plot:
                    xmin = 1e-3
                    fig, ax = plt.subplots()
                    
                    xpl = np.linspace(xmin, 1-xmin, 1000)
                    ax.plot(xpl, xpl, lw=1.5, color='k')

                    ax.scatter(af1, af2)
                    ax.grid(True)
                    ax.set_xscale('logit')
                    ax.set_yscale('logit')
                    ax.set_xlim(xmin, 1-xmin)
                    ax.set_ylim(xmin, 1-xmin)
                    ax.set_xlabel('Earlier sample '+samplename1)
                    ax.set_ylabel('Later sample '+samplename2)
                    ax.set_title('Interval: '+str(dt)+' days')

                    plt.tight_layout()

                    plt.ion()
                    plt.show()


    if use_plot:

        fig, ax = plt.subplots()
        keys = data['eta'].keys()
        y, x, ifrs = zip(*[(data['eta'][key], data['dt'][key], int(key[1][1])-1)
                         for key in keys])
        color = cm.jet(1.0 * np.array(ifrs) / 6)

        xpl = np.logspace(0, 3, 100)
        h = ax.plot(xpl, 2e-3 + 2e-4 * xpl**1, lw=1.5, c='k')[0]
        ax.scatter(x, y, color=color)

        ax.set_xlabel('Time interval [days]')
        ax.set_ylabel('Error magnitude: var / x (1 - x)')
        
        ax.grid(True)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(1, 2000)
        ax.set_ylim(1e-3, 1)

        plt.tight_layout()

        plt.ion()
        plt.show()

    sys.exit()
