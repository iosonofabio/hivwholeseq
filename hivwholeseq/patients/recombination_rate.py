# vim: fdm=marker
'''
author:     Fabio Zanini
date:       20/03/14
content:    Measure recombination in various ways.
'''
# Modules
import sys
import os
import argparse
from itertools import izip
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.samples import SamplePat



# Functions
def broadcast_one_site_to_two_sites(array):
    '''Broadcast a one-site statistics for two-sites'''
    array_bc = np.tile(array, np.prod(array.shape))
    array_bc = array_bc.reshape((array.shape[0], array.shape[1],
                                 array.shape[0], array.shape[1]))
    array_bc = array_bc.swapaxes(1, 2)
    return array_bc



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Filter mapped reads')
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot results')

    args = parser.parse_args()
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    use_plot = args.plot


    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]

    ds_all = []
    ys_all = []
    for fragment in fragments:
        for pname, patient in patients.iterrows():
            if VERBOSE >= 1:
                print pname, fragment

            patient = Patient(patient)

            refseq = patient.get_reference(fragment)

            if VERBOSE >= 2:
                print 'Get allele frequency trajectories'
            aft, indt = patient.get_allele_frequency_trajectories(fragment)
            times = patient.times[indt]

            if VERBOSE >= 2:
                print 'Get initial cocounts'
            sample0 = SamplePat(patient.samples.iloc[0])
            # FIXME: complete those matrices
            try:
                acc = sample0.get_allele_cocounts(fragment)
            except IOError:
                ds_all.append([])
                ys_all.append([])
                continue

            # FIXME: Einstein's style sums should be equivalent to tile + swapaxes
            # EXCEPT that einsum ignores the mask I should DOUBLE CHECK that. Btw,
            # if speed is an issue and einsum is faster, we can just do it twice,
            # with data and mask, and then reglue those togetherfaster, we can just do it twice,
            # with data and mask, and then reglue those together..

            if VERBOSE >= 2:
                print 'Restrict to polymorphic sites'
            # FIXME: deal better with low viral load
            ind_poly = (aft > 0.01).any(axis=0) & (aft < 1 - 0.01).any(axis=0)
            ind_poly_cc = broadcast_one_site_to_two_sites(ind_poly)
            ind_poly_cc = ind_poly_cc & ind_poly_cc.swapaxes(0, 1).swapaxes(2, 3)

            if VERBOSE >= 2:
                print 'Restrict to pairs with enough cocoverage'
            # FIXME: make a variable threshold
            ind_cocoverage = (acc.sum(axis=0).sum(axis=0) > 100)
            ind_cocoverage_cc = np.einsum('ij,kl->klij',
                                          ind_cocoverage,
                                          np.ones(acc.shape[:2], bool))

            if VERBOSE >= 2:
                print 'Calculate coallele frequencies from cocounts, excluding low-coverage and low-frequency'
            acf = 1.0 * acc / acc.sum(axis=0).sum(axis=0)
            acf[ind_cocoverage_cc & (acf < 1e-2)] = 0
            # Renormalize (small effect)
            acf /= acf.sum(axis=0).sum(axis=0)

            if VERBOSE >= 2:
                print 'Take a subset of pairs'
            indi_good_cc = np.transpose((ind_poly_cc & ind_cocoverage_cc).nonzero())
            indi_ss = indi_good_cc[np.random.randint(len(indi_good_cc), size=1000)]
            indi_ts = np.random.randint(len(times), size=(len(indi_ss), 2))

            if VERBOSE >= 2:
                print 'Computing statistics'
            ds = np.zeros(len(indi_ss))
            ys = np.zeros_like(ds)
            for i, (indi, it) in enumerate(izip(indi_ss, indi_ts)):
                it1 = np.min(it)
                it2 = np.max(it)

                if it1 == it2:
                    continue

                # NOTE: tuples are special for indexing (!)
                sgn = np.sign(acf[tuple(indi)] - \
                              aft[0, indi[0], indi[2]] * aft[0, indi[1], indi[3]])
                daf1 = aft[it2, indi[0], indi[2]] - aft[it1, indi[0], indi[2]]
                daf2 = aft[it2, indi[1], indi[3]] - aft[it1, indi[1], indi[3]]
                y = sgn * daf1 * daf2

                dt = times[it2] - times[it1]
                dpos = np.abs(indi[3] - indi[2])
                d = dpos * dt

                if VERBOSE >= 4:
                    print acc[:, :, indi[2], indi[3]]
                    print acf[:, :, indi[2], indi[3]]
                    import ipdb; ipdb.set_trace()

                ds[i] = d
                ys[i] = y

            ds_all.append(ds)
            ys_all.append(ys)

    if use_plot:

        if VERBOSE >= 2:
            print 'Plotting'

        dsa = np.concatenate(ds_all)
        ysa = np.concatenate(ys_all)

        bins = np.linspace(0.001, dsa.max() + 0.001, 30)
        x = 0.5 * (bins[1:] + bins[:-1])
        y = np.array([ysa[(dsa >= bins[i]) & (dsa < bins[i+1])].mean()
                      for i in xrange(len(bins) - 1)])
        dy = np.array([ysa[(dsa >= bins[i]) & (dsa < bins[i+1])].std()
                       for i in xrange(len(bins) - 1)])

        fig, ax = plt.subplots()
        ax.set_xlabel('Distance [days * bp]')
        ax.set_ylabel('sgn$C_0 \cdot d\\nu_1 \cdot d\\nu_2$')
        ax.grid(True)

        ax.errorbar(x, y, dy, lw=2)
        
        plt.tight_layout()
        plt.ion()
        plt.show()







        #if VERBOSE >= 1:
        #    print 'Initializing matrix of coallele frequencies'
        #coaft = np.empty((len(samplenames), len(alpha), len(alpha), len(refseq), len(refseq)), float)
        #if VERBOSE >= 1:
        #    print 'Collecting cocounts and normalizing'
        #for i, samplename in enumerate(samplenames):
        #    if VERBOSE >= 2:
        #        print pname, fragment, samplename
        #    cocounts = np.load(get_allele_cocounts_filename(pname, samplename, fragment))
        #    # NOTE: uncovered will get 0 correlation
        #    coaft[i] = 1.0 * cocounts / (cocounts.sum(axis=0).sum(axis=0) + 0.1)
        #    del cocounts

        #if VERBOSE >= 1:
        #    print 'Getting allele frequencies'
        #aft = np.load(get_allele_frequency_trajectories_filename(pname, fragment))[:coaft.shape[0]]
        #if VERBOSE >= 1:
        #    print 'Broadcasting allele frequency product (outer tensor product)'
        #aftbroad = np.einsum('mij...,m...kl->mikjl', aft, np.ones_like(aft))
        #aftbroadT = aftbroad.swapaxes(1, 2).swapaxes(3, 4)
        #if VERBOSE >= 1:
        #    print 'Subtract product of allele frequencies'
        #LD = coaft - (aftbroad * aftbroadT)

        ## Mask LD when allele freqs are too close to 0 or 1
        #LD = np.ma.array(LD)
        #LD[(aftbroad < 1e-2) | (aftbroad > 1 - 1e-2) | \
        #   (aftbroadT < 1e-2) | (aftbroadT > 1 - 1e-2)] = np.ma.masked

        ## Take the correlation coefficient, i.e. divide by sqrt(pi qi pj qj)
        #corr = LD / np.sqrt(aftbroad * (1 - aftbroad) * aftbroadT * (1 - aftbroadT))

        #if VERBOSE >= 1:
        #    print 'Finding initial highly correlated pairs'

        ## FIXME: only use this criterion because we have t0 twice
        #ind_highLD0 = (corr[0] > 1).nonzero()
        ##ind_highLD0 = ((LD[0] > 0.1) & (LD[1] > 0.1) & (aftbroad[0] > 1e-3) & (aftbroad[0] < 0.5)).nonzero()
        #pairs_high_LD0 = zip(*ind_highLD0)
        #LD_pairs = np.array([LD[:, pair[0], pair[1], pair[2], pair[3]] for pair in pairs_high_LD0]).T
        #corr_pairs = np.array([corr[:, pair[0], pair[1], pair[2], pair[3]] for pair in pairs_high_LD0]).T
        #print corr_pairs.T
