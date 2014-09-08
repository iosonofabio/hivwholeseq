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
from hivwholeseq.patients.patients import load_patient, convert_date_deltas_to_float
from hivwholeseq.patients.filenames import get_allele_frequency_trajectories_filename, \
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
    parser = argparse.ArgumentParser(description='Get divergence and diversity',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the allele frequency trajectories')
    parser.add_argument('--PCR1', type=int, default=1,
                        help='Take only PCR1 samples [0=off, 1=where both available, 2=always]')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose
    plot = args.plot
    use_PCR1 = args.PCR1

    patient = load_patient(pname)
    patient.discard_nonsequenced_samples()

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    
    if plot:
        fig, ax = plt.subplots(1, 1)
        ax.set_xlabel('Time from transmission [days]')
        ax.set_ylabel('Divergence/diversity')
        ax.set_yscale('log')
        ax.set_title(pname)

    dgs = {}
    dss = {}
    for ifr, fragment in enumerate(fragments):
        if VERBOSE >= 1:
            print pname, fragment

        aft, ind = patient.get_allele_frequency_trajectories(fragment, use_PCR1=use_PCR1)
        times = patient.times[ind]

        cons_ind = aft[0].argmax(axis=0)
        dg = np.zeros(len(times))
        ds = np.zeros_like(dg)
        for pos in xrange(aft.shape[2]):
            af = aft[:, :, pos]
            af_nonanc = 1.0 - af[:, cons_ind[pos]]
            dg += af_nonanc
            ds += (af * (1 - af)).sum(axis=1)

        dg /= aft.shape[2]
        ds /= aft.shape[2]

        dgs[fragment] = dg
        dss[fragment] = ds

        if plot:
            ax.plot(times, dg, lw=2, label=fragment,
                    color=cm.jet(1.0 * ifr / len(fragments)))

            ax.plot(times, ds, lw=2, ls='--',
                    color=cm.jet(1.0 * ifr / len(fragments)))

    if plot:
        ax.legend(loc=4, fontsize=12)

    plt.ion()
    plt.show()

