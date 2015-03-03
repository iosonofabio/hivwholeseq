# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/09/14
content:    Get divergence and diversity of the patient.
'''
# Modules
import os
import argparse
from operator import itemgetter
import numpy as np
from Bio import SeqIO
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.patients import load_patients, Patient



# Functions
def get_divergence(aft):
    '''Get divergence from allele frequency trajectories'''
    cons_ind = Patient.get_initial_consensus_noinsertions(aft, return_ind=True)
    dg = 1 - aft[:, cons_ind, np.arange(aft.shape[2])].mean(axis=1)
    return dg


def get_diversity(aft):
    '''Get diversity from allele frequency trajectories'''
    ds = (aft * (1 - aft)).sum(axis=1).mean(axis=1)
    return ds



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get divergence and diversity',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Regions to analyze (e.g. F1 V3)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the allele frequency trajectories')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    plot = args.plot

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    data = []

    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()

        for ifr, region in enumerate(regions):
            if VERBOSE >= 1:
                print pname, region

            aft, ind = patient.get_allele_frequency_trajectories(region,
                                                                 cov_min=10)
            times = patient.times[ind]

            dg = get_divergence(aft)
            ds = get_diversity(aft)

            data.append({'pname': pname, 'region': region, 'dg': dg, 'ds': ds, 't': times})

    if plot:
        fig, ax = plt.subplots(1, 1)
        ax.set_xlabel('Time from transmission [days]')
        ax.set_ylabel('Divergence [solid]\nDiversity [dashed]')
        #ax.set_yscale('log')

        for i, d in enumerate(data):
            pname = d['pname']
            region = d['region']
            dg = d['dg']
            ds = d['ds']
            times = d['t']

            ax.plot(times, dg, lw=2, label=', '.join([pname, region]),
                    color=cm.jet(1.0 * i / len(data)))

            ax.plot(times, ds, lw=2, ls='--',
                    color=cm.jet(1.0 * i / len(data)))

        ax.legend(loc=4, fontsize=12)
        ax.grid(True)
        plt.tight_layout()

        plt.ion()
        plt.show()

