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
from hivwholeseq.patients.one_site_statistics import \
        plot_allele_frequency_trajectories_from_counts as plot_nus_from_act
from hivwholeseq.patients.one_site_statistics import \
        plot_allele_frequency_trajectories_from_counts_3d as plot_nus_from_act_3d
from hivwholeseq.fork_cluster import fork_get_allele_frequency_trajectory as fork_self


# Functions
def get_divergence(aft):
    '''Get divergence from allele frequency trajectories'''
    cons_ind = Patient.get_initial_consensus_noinsertions(aft)
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
    parser.add_argument('--patients', nargs='+', default=['all'],
                        help='Patients to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the allele frequency trajectories')
    parser.add_argument('--saveplot', action='store_true',
                        help='Save the plot to file')
    parser.add_argument('--PCR1', type=int, default=1,
                        help='Take only PCR1 samples [0=off, 1=where both available, 2=always]')

    args = parser.parse_args()
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    plot = args.plot
    saveplot = args.saveplot
    use_PCR1 = args.PCR1

    patients = load_patients()
    if pnames != ['all']:
        patients = patients.iloc[patients.index.isin(pnames)]

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    dgs = {}
    dss = {}
    
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()

        if plot:
            fig, ax = plt.subplots(1, 1)
            ax.set_xlabel('Time from transmission [days]')
            ax.set_ylabel('Divergence [solid]\nDiversity [dashed]')
            ax.set_yscale('log')
            ax.set_title(pname)

        for ifr, fragment in enumerate(fragments):
            if VERBOSE >= 1:
                print pname, fragment

            aft, ind = patient.get_allele_frequency_trajectories(fragment,
                                                                 use_PCR1=use_PCR1,
                                                                 cov_min=10)
            times = patient.times[ind]

            dg = get_divergence(aft)
            ds = get_diversity(aft)

            dgs[(pname, fragment)] = dg
            dss[(pname, fragment)] = ds

            if plot:
                ax.plot(times, dg, lw=2, label=fragment,
                        color=cm.jet(1.0 * ifr / len(fragments)))

                ax.plot(times, ds, lw=2, ls='--',
                        color=cm.jet(1.0 * ifr / len(fragments)))

        if plot:
            ax.legend(loc=4, fontsize=12)
            ax.grid(True)
            plt.tight_layout()

            if saveplot:
                plt.savefig('/ebio/ag-neher/home/fzanini/phd/sequencing/figures/'+\
                            'divergence_diversity_'+pname+'.png')
                plt.close(fig)

    if not saveplot:
        plt.ion()
        plt.show()

