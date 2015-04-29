# vim: fdm=marker
'''
author:     Fabio Zanini
date:       02/07/14
content:    Explore sites that are very conserved within our patient set, across
            the whole infection, compared to their behaviour within the subtype.
'''
# Modules
import os
import sys
import argparse
from itertools import izip, combinations, chain
from collections import defaultdict
from operator import itemgetter
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt
from Bio import SeqIO, AlignIO
from scipy.stats import pearsonr, spearmanr

from hivwholeseq.miseq import alpha, alphal
from hivwholeseq.patients.filenames import root_patient_folder
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.cross_sectional.get_subtype_consensus import (
    get_subtype_reference_alignment_consensus)
from hivwholeseq.cross_sectional.get_subtype_entropy import (
    get_subtype_reference_alignment_entropy)
from hivwholeseq.cross_sectional.get_subtype_allele_frequencies import (
    get_subtype_reference_alignment_allele_frequencies)



# Globals
colors = {'p17': 'r',
          'p24': 'g',
          'PR': 'r',
          'RT': 'g',
          'p15': 'purple',
          'IN': 'orange',
          'gp41': 'r'}
regions = ['p17', 'p24', 'PR', 'RT', 'p15', 'IN', 'gp41']
pnames = ['20097', '15363', '15823', '9669', '20529', '15241', '15376', '15319']




# Functions
def get_allele_freqs_alignment(alim, positions=None, VERBOSE=0):
    '''Get allele frequencies from alignment'''
    from hivwholeseq.miseq import alpha, alphal

    # Make sure it's a numpy matrix, efficiently
    alim = np.asarray(alim, 'S1')

    if positions is None:
        positions = np.arange(alim.shape[-1])

    # Get allele freqs ignoring ambiguous: N, Y, R, etc. and renormalizing
    num = np.array([(alim == a).mean(axis=0) for a in alpha[:5]])
    num /= num.sum(axis=0)

    return num[:, positions]


def get_entropy_pats(afts, VERBOSE=0):
    '''Calculate entropies'''
    n_patients = afts.shape[0]
    lseq = afts.shape[1]

    if VERBOSE >= 1:
        print 'Calculate entropy for each patient'
    S = -np.ones((n_patients, lseq), float)
    S_time = np.ma.zeros((n_patients, lseq), object)
    times_covs = np.zeros(n_patients, object)
    for k, aft_pat in enumerate(afts):
        times_cov = aft_pat[0, 0] >= -0.1
        times_covs[k] = times_cov
        if times_cov.sum() < 2:
            S_time[k, :] = np.ma.masked
            continue
        
        for pos, aft_pos in enumerate(aft_pat):
            # FIXME: deal better with missing data (low coverage)
            Stmp = 0
            Stmp_time = np.zeros(len(times_cov))
            for j, aft_nuc in enumerate(aft_pos):
                aft_nuc_cov = aft_nuc[times_cov]
                Stmp -= ((aft_nuc_cov + 1e-8) * np.log2(aft_nuc_cov + 1e-8)).mean()
                Stmp_time -= ((aft_nuc + 1e-8) * np.log2(aft_nuc + 1e-8))
            S[k, pos] = Stmp
            S_time[k, pos] = Stmp_time

    return {'mean': S,
            'time': S_time}


def collect_data(pnames, regions, VERBOSE=0):
    '''Collect data'''

    data = []
    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    for region in regions:
        if VERBOSE >= 1:
            print region

        if VERBOSE >= 2:
            print 'Get subtype consensus'
        conssub = get_subtype_reference_alignment_consensus(region, VERBOSE=VERBOSE)

        if VERBOSE >= 2:
            print 'Get subtype allele frequencies'
        af_sub = get_subtype_reference_alignment_allele_frequencies(region)

        if VERBOSE >= 2:
            print 'Get subtype entropy'
        Ssub = get_subtype_reference_alignment_entropy(region, VERBOSE=VERBOSE)

        for ipat, (pname, patient) in enumerate(patients.iterrows()):
            pcode = patient.code
            if VERBOSE >= 2:
                print pname, pcode

            patient = Patient(patient)
            aft, ind = patient.get_allele_frequency_trajectories(region,
                                                                 cov_min=1000,
                                                                 depth_min=300,
                                                                 VERBOSE=VERBOSE)
            if len(ind) == 0:
                if VERBOSE >= 2:
                    print 'No time points found: skip.'
                continue

            times = patient.times[ind]

            if VERBOSE >= 2:
                print 'Get coordinate map'
            coomap = patient.get_map_coordinates_reference(region, refname=('HXB2', region))

            if VERBOSE >= 2:
                print 'Get initial consensus'
            icons = patient.get_initial_consensus_noinsertions(aft, VERBOSE=VERBOSE,
                                                               return_ind=True)
            consm = alpha[icons]

            # Get the map as a dictionary from patient to subtype
            coomapd = dict(coomap[:, ::-1])

            if VERBOSE >= 2:
                print 'Get alleles'
            for posdna in xrange(aft.shape[-1]):
                if VERBOSE >= 3:
                    print posdna

                # Look for this position in the subtype alignment
                if posdna not in coomapd:
                    continue
                pos_sub = coomapd[posdna]
                if (pos_sub >= af_sub.shape[1]):
                    continue

                # Ancestral allele
                ianc = icons[posdna]
                anc = alpha[ianc]

                # Discard if the initial time point is already polymorphic
                aft_anc0 = aft[0, ianc, posdna]
                if aft_anc0 < 0.9:
                    continue


                # Get only non-masked time points
                aft_der = aft[:, 0, posdna]
                indpost = -aft_der.mask
                if indpost.sum() == 0:
                    continue
                timespos = times[indpost]
                aft_nm = aft[indpost, :, posdna]

                Spos_sub = Ssub[pos_sub]
                conspos_sub = conssub[pos_sub]

                if (anc == conspos_sub):
                    away_conssub = 'away'
                else:
                    away_conssub = 'to'

                for it, time in enumerate(timespos):
                    S = -(aft_nm[it] * np.log2(aft_nm[it] + 1e-10)).sum()
                    data.append((region, pcode,
                                 anc,
                                 Spos_sub, conspos_sub, away_conssub,
                                 time, S))

    data = pd.DataFrame(data=data, columns=('region', 'pcode',
                                            'anc',
                                            'Ssub', 'conssub', 'awayto',
                                            'time', 'S'))


    return data


def calculate_correlation(data, VERBOSE=0):
    '''Calculate correlation coefficients'''
    from scipy.stats import pearsonr, spearmanr

    datap = []
    for (pcode, region), datum in data.groupby(['pcode', 'region']):
        for t, datumt in datum.groupby('time'):
            M = np.array(datumt[['Ssub', 'S']]).T
            rho, P = spearmanr(M[0], M[1])
            datap.append({'pcode': pcode,
                          'region': region,
                          'time': t,
                          'rho': rho,
                         })
    datap = pd.DataFrame(datap)

    return datap


def plot_correlation(datap, VERBOSE=0, colormap='jet'):
    '''Plot the correlation between intrapatient and subtype diversity'''
    if isinstance(colormap, basestring):
        if colormap == 'jet':
            colormap = cm.jet
        else:
            raise ValueError('colormap not recognized')

    pcodes = datap['pcode'].unique().tolist()
    for region, datump in datap.groupby('region'):

        fig, ax = plt.subplots()


        for ip, (pcode, datumpp) in enumerate(datump.groupby('pcode')):
            x = np.array(datumpp['time'])
            y = np.array(datumpp['rho'])
            ax.plot(x, y,
                    lw=2,
                    color=colormap(1.0 * ip / len(pcodes)),
                   )

        ax.set_xlabel('Time [days from infection')
        ax.set_ylabel('Spearmanr pat/subtype B')
        ax.grid(True)



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Explore conservation levels across patients and subtype',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+', default=pnames,
                        help='Patient to analyze')
    parser.add_argument('--regions', nargs='+', default=regions,
                        help='Regions to analyze (e.g. F1 p17)')
    parser.add_argument('--verbose', type=int, default=2,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot results')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    plot = args.plot

    data = collect_data(pnames, regions, VERBOSE=VERBOSE)

    datap = calculate_correlation(data, VERBOSE=VERBOSE)


    if plot:
        plot_correlation(datap, VERBOSE=VERBOSE)
        plt.ion()
        plt.show()


