# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/09/14
content:    Get divergence and diversity of the patient.
'''
# Modules
import os
import sys
import argparse
from operator import itemgetter
import numpy as np
from Bio import SeqIO
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha
from hivwholeseq.utils.argparse import PatientsAction
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.filenames import get_allele_frequency_trajectories_filename, \
        get_allele_count_trajectories_filename
from hivwholeseq.patients.one_site_statistics import \
        plot_allele_frequency_trajectories_from_counts as plot_nus_from_act
from hivwholeseq.patients.one_site_statistics import \
        plot_allele_frequency_trajectories_from_counts_3d as plot_nus_from_act_3d
from hivwholeseq.patients.one_site_statistics import get_allele_count_trajectories
from hivwholeseq.cluster.fork_cluster import fork_get_allele_frequency_trajectory as fork_self


# Functions
def get_divergence_trajectory_local(pname, fragment, block_length=150, VERBOSE=0):
    '''Get local divergence trajectory'''
    from hivwholeseq.patients.filenames import get_divergence_trajectories_local_filename
    fn = get_divergence_trajectories_local_filename(pname, fragment)
    npz = np.load(fn)
    return (npz['dg'], npz['ind'], npz['block_length'], npz['L'])


def get_diversity_trajectory_local(pname, fragment, block_length=150, VERBOSE=0):
    '''Get local diversity trajectory'''
    from hivwholeseq.patients.filenames import get_diversity_trajectories_local_filename
    fn = get_diversity_trajectories_local_filename(pname, fragment)
    npz = np.load(fn)
    return (npz['ds'], npz['ind'], npz['block_length'], npz['L'])


def plot_divdiv_trajectory(patient, VERBOSE=0):
    '''Plot the trajectory of divergence and diversity'''
    ind = patient.ind
    dg = patient.dg
    ds = patient.ds
    L = patient.L
    block_length = patient.block_length

    x = (0.5 + np.arange(dg.shape[1])) * block_length
    times = patient.times[ind]

    if hasattr(patient, 'covt'):
        covt = patient.covt
        use_coverage = True
    else:
        use_coverage = False

    fig, axs = plt.subplots(1, 2 + use_coverage,
                            figsize=(15 + 4 * use_coverage, 8))
    fig.suptitle(pname+', '+fragment+\
                 ' (block length = '+str(block_length)+')', fontsize=18)
    ax1 = axs[0]
    ax1.set_xlabel('Position [bp]')
    ax1.set_ylabel('Divergence')
    ax1.set_yscale('log')
    ax1.set_ylim(0.9e-4, 5e-1)
    ax1.set_xlim(-150, L + 50)
    ax1.grid(True)

    ax2 = axs[1]
    ax2.set_xlabel('Position [bp]')
    ax2.set_ylabel('Diversity')
    ax2.set_yscale('log')
    ax2.set_ylim(0.9e-4, 5e-1)
    ax2.set_xlim(-150, L + 50)
    ax2.grid(True)

    if use_coverage:
        ax3 = axs[2]
        ax3.set_xlabel('Position [bp]')
        ax3.set_ylabel('Coverage')
        ax3.set_yscale('log')
        ax3.set_ylim(1e-1, 7e5)
        ax3.set_xlim(-150, L + 50)
        ax3.grid(True)

    for it in xrange(len(times)):
        ax1.plot(x, dg[it] + 1e-4,
                 lw=2, label=str(times[it]),
                 color=cm.jet(1.0 * it / len(times)),
                 alpha=0.5)

        ax2.plot(x, ds[it] + 1e-4,
                 lw=2,
                 color=cm.jet(1.0 * it / len(times)),
                 alpha=0.5)

        if use_coverage:
            ax3.plot(np.arange(L), covt[it],
                     lw=2, label=patient.samples[ind].index[it],
                     color=cm.jet(1.0 * it / len(times)))

    # Plot averages
    ax1.scatter([-100] * len(times), dg.mean(axis=1), s=40,
                c=[cm.jet(1.0 * it / len(times))
                   for it in xrange(len(times))],
                edgecolor='none', alpha=0.5)

    ax2.scatter([-100] * len(times), ds.mean(axis=1), s=40,
                            c=[cm.jet(1.0 * it / len(times))
                               for it in xrange(len(times))],
                            edgecolor='none', alpha=0.5)

    if use_coverage:
        ax3.scatter([-100] * len(times), covt.mean(axis=1), s=40,
                    c=[cm.jet(1.0 * it / len(times)) for it in xrange(len(times))],
                    edgecolor='none', alpha=0.5)

    ax1.legend(loc=4, fontsize=10, ncol=(1 + ((len(times) - 1) // 3)),
               title='Time from transmission [days]')

    if use_coverage:
        ax3.legend(loc=3, fontsize=10, ncol=(1 + ((len(times) - 1) // 3)),
                   title='Sample name:')

    if fragment == 'genomewide':
        refann = patient.get_reference(fragment, 'gb')
        for fea in refann.features:
            if fea.type != 'fragment':
                continue

            frn = int(fea.id[1])
            frn_even = 1 - (frn % 2)
            x = (fea.location.nofuzzy_start, fea.location.nofuzzy_end)
            y = [0.95**(1 + frn_even) * ax1.get_ylim()[1]] * 2
            ax1.plot(x, y, color='k', lw=2)

    plt.tight_layout(rect=(0, 0, 1, 0.95))



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get local divergence and diversity',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', action=PatientsAction,
                        help='Patients to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the result')
    parser.add_argument('--block-length', type=int, default=150,
                        help='Length of block to consider')
    parser.add_argument('--include-coverage', action='store_true', dest='include_cov',
                        help='Also calculate/plot the coverage')

    args = parser.parse_args()
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    plot = args.plot
    block_length = args.block_length
    use_coverage = args.include_cov

    patients = load_patients()
    if pnames is not None:
        patients = patients.iloc[patients.index.isin(pnames)]

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()

        for ifr, fragment in enumerate(fragments):
            if VERBOSE >= 1:
                print pname, fragment

            dg, ind, block_length, L = \
                    patient.get_divergence_trajectory_local(fragment,
                                                            block_length=block_length)
            ds, ind, block_length, L = \
                    patient.get_diversity_trajectory_local(fragment,
                                                           block_length=block_length)
            patient.dg = dg
            patient.ds = ds
            patient.ind = ind
            patient.L = L
            patient.block_length = block_length

            if use_coverage:
                (covt, ind2) = patient.get_coverage_trajectories(fragment)
                if set(ind).symmetric_difference(set(ind2)):
                    raise ValueError('Indices for allele freqs and coverage differ!')
                patient.covt = covt

            if plot:
                plot_divdiv_trajectory(patient, VERBOSE=VERBOSE)

       
    if plot:   
        plt.ion()
        plt.show()

