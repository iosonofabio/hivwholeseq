# vim: fdm=marker
'''
author:     Fabio Zanini
date:       19/05/14
content:    Calculate and plot the propagator of allele frequencies.
'''
# Modules
import os
import sys
import argparse
from operator import itemgetter
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha
from hivwholeseq.reference import load_custom_reference, load_custom_alignment
from hivwholeseq.patients.patients import load_patients, filter_patients_n_times, Patient
from hivwholeseq.patients.filenames import get_initial_reference_filename
from hivwholeseq.patients.one_site_statistics import get_allele_count_trajectories
from hivwholeseq.utils.import plot
from hivwholeseq.patients.propagator_allele_frequency import Propagator, \
        plot_propagator_theory
from hivwholeseq.patients.plot_SFS_by_entropy import get_coordinate_map



# Script
if __name__ == '__main__': 

    parser = argparse.ArgumentParser(description='Propagator for allele frequencies',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save the propagator to file')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the propagator')
    parser.add_argument('--deltat', type=int, nargs=2, default=[100, 300],
                        help='Time in days between final and initial (range)')
    parser.add_argument('--logit', action='store_true',
                        help='use logit scale (log(x/(1-x)) in the plots')
    parser.add_argument('--min-depth', type=int, default=100, dest='min_depth',
                        help='Minimal depth to consider the site')

    args = parser.parse_args()
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    save_to_file = args.save
    plot = args.plot
    dt = args.deltat
    use_logit = args.logit
    depth_min = args.min_depth

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    # Select type M entropy strata
    S_bins = np.array([0.005, 0.05, 0.5, 2])

    if VERBOSE >= 1:
        print 'Load alignment, reference, and coordinate map'
    #ali = load_custom_alignment('HIV1_FLT_2013_genome_DNA')
    #alim = np.array(ali, 'S1')
    S = np.zeros(alim.shape[1])
    for a in alpha[:5]:
        nu = (alim == a).mean(axis=0)
        S -= nu * np.log2(nu + 1e-8)

    refname = 'HXB2'
    refseq = load_custom_reference('HXB2', format='gb')
    mapali = get_coordinate_map(refname, ali)
    if len(refseq) != mapali[0, -1] + 1:
        raise ValueError('Reference '+refname+' in alignment is not complete')
    Sref = S[mapali[1]]

    Srefind = - np.ones(len(Sref), int)
    for i, Sbin in enumerate(S_bins[:-2]):
        ind = (Sref >= Sbin) & (Sref < S_bins[i + 1])
        Srefind[ind] = i
    Srefind[Srefind < 0] = len(S_bins) - 2

    # Prepare output structures
    n_binsx = 5
    binsy = [0.,
             0.002,
             0.01,
             0.05,
             0.2,
             0.5,
             0.8,
             0.95,
             0.99,
             0.998,
             1.]
    pps = [Propagator(n_binsx, binsy=binsy, use_logit=use_logit) for _ in S_bins[:-1]]

    histrs = [pp.histogram for pp in pps]

    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        samplenames = patient.samples.index

        if not fragments:
            fragments = ['F'+str(i) for i in xrange(1, 7)]
        if VERBOSE >= 2:
            print 'fragments', fragments
    
        # Iterate over samples and fragments
        for fragment in fragments:
            if VERBOSE >= 1:
                print pname, fragment
    
            mapco = patient.get_map_coordinates_reference(fragment, refname=refname)
    
            aft, ind = patient.get_allele_frequency_trajectories(fragment,
                                                                 cov_min=depth_min)

            ts = patient.times[ind]
    
            # Collect counts
            for i in xrange(aft.shape[0] - 1):
                for j in xrange(i + 1, aft.shape[0]):
                    dij = j - i
                    if ((ts[j] - ts[i]) > dt[1]) or ((ts[j] - ts[i]) < dt[0]):
                        continue

                    for (pos_ref, pos) in mapco:
                        k = Srefind[pos_ref]
                        histrs[k] += np.histogram2d(aft[i, :, pos],
                                                    aft[j, :, pos],
                                                    bins=[pps[k].binsx, pps[k].binsy])[0]

    if plot:
        for i, pp in enumerate(pps):
            if VERBOSE >= 1:
                print 'Plot HIV propagator: S e ['+str(S_bins[i])+', '+str(S_bins[i+1])+']'
            title = 'Propagator for allele frequencies\n'+\
                    '$\Delta t = '+str(dt)+'$ days, '+str(fragments)+'\n'+\
                    '$S \in ['+str(S_bins[i])+', '+str(S_bins[i+1])+']$'
            pp.plot(title=title, heatmap=False)

        if VERBOSE >= 1:
            print 'Calculate and plot theory'

        t = 1.0 * np.mean(dt) / 500
        xis = pp.binsxc
        plot_propagator_theory(xis, t, model='BSC',
                               logit=use_logit, xlim=[pp.binsyc[0], pp.binsyc[-1]])
        plot_propagator_theory(xis, t, model='neutral',
                               logit=use_logit, xlim=[pp.binsyc[0], pp.binsyc[-1]])
        
        plt.ion()
        plt.show()
