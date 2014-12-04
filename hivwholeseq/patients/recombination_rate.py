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


    # data contains lists of:
    # distance, dt, dbp, dymax, dymin, sgn(cov0), afmax, ipat, obs
    idata = 0 
    data = np.zeros((2000000, 9), float)
    for ipat, (pname, patient) in enumerate(patients.iterrows()):
        patient = Patient(patient)

        for fragment in fragments:
            if VERBOSE >= 1:
                print pname, fragment

            refseq = patient.get_reference(fragment)

            if VERBOSE >= 2:
                print 'Get allele frequency trajectories'
            aft, indt = patient.get_allele_frequency_trajectories(fragment)
            times = patient.times[indt]

            if VERBOSE >= 2:
                print 'Iterate over initial times'

            for it0, t0 in enumerate(times[:-1]):
                if VERBOSE >= 2:
                    print 'Get initial cocounts for time point', it0 + 1, 'out of', len(times)
                sample0 = SamplePat(patient.samples.iloc[it0])
                try:
                    acc = sample0.get_allele_cocounts(fragment)
                except IOError:
                    continue

                if VERBOSE >= 2:
                    print 'Restrict to polymorphic sites and minor alleles:',
                # FIXME: deal better with low viral load
                af0 = aft[it0]
                ind_poly = (af0 > 0.05) & (af0 < 0.5)
                indi_poly = ind_poly.nonzero()
                n_poly = len(indi_poly[0])
                if VERBOSE >= 2:
                    print n_poly, 'found'

                if VERBOSE >= 2:
                    print 'Iterating over pairs'
                n_pairs_checked = 0
                for ip1 in xrange(n_poly - 1):
                    for ip2 in xrange(ip1 + 1, n_poly):
                        if VERBOSE >= 3:
                            print 'Pairs checked:', n_pairs_checked
                        n_pairs_checked += 1

                        pos1, pos2 = indpos12 = tuple(indi_poly[1][[ip1, ip2]])
                        if pos1 == pos2:
                            continue

                        indp1 = (indi_poly[0][ip1], indi_poly[1][ip1])
                        indp2 = (indi_poly[0][ip2], indi_poly[1][ip2])
                        indp12 = (indi_poly[0][ip1], indi_poly[0][ip2],
                                  indi_poly[1][ip1], indi_poly[1][ip2])

                        if VERBOSE >= 3:
                            print 'Check for pair cocoverage:',
                        accpos12 = acc[:, :, indpos12[0], indpos12[1]]
                        cocov = accpos12.sum()
                        if VERBOSE >= 3:
                            print cocov
                        if cocov < 100:
                            continue

                        if VERBOSE >= 3:
                            print 'Calculate coallele frequencies:',
                        acf12 = 1.0 * acc[indp12] / cocov
                        if VERBOSE >= 3:
                            print acf12

                        if VERBOSE >= 3:
                            print 'Check initial correlation:',
                        corr12 = (acf12 - af0[indp1] * af0[indp2]) / (af0[indp1] * af0[indp2])
                        sgncorr12 = np.sign(corr12)
                        if VERBOSE >= 3:
                            print corr12
                            print 'The sign of the correlation is', sgncorr12

                        if VERBOSE >= 3:
                            print 'Take next time point'
                        it1 = it0 + 1
                        t1 = times[it1]
                        af1 = aft[it1]

                        if VERBOSE >= 3:
                            print 'Calculate delta allele freqs'
                        daf1 = af1[indp1] - af0[indp1]
                        daf2 = af1[indp2] - af0[indp2]

                        if VERBOSE >= 3:
                            print 'Calculate spatiotemporal distance:',
                        d = np.abs((t1 - t0) * (pos2 - pos1))
                        if VERBOSE >= 3:
                            print d, 'days * bp'

                        if VERBOSE >= 3:
                            print 'Calculate linkage statistics:',
                        y = sgncorr12 * daf1 * daf2
                        if VERBOSE >= 3:
                            print y

                        # The first af if the largest changing
                        dafs = (daf1, daf2)
                        daf_max_ind = np.argmax(map(np.abs, dafs))
                        daf_min_ind = int(not bool(daf_max_ind))
                        daf_max = dafs[daf_max_ind]
                        daf_min = dafs[daf_min_ind]
                        data[idata] = (d, t1 - t0, np.abs(pos2 - pos1),
                                       daf_max, daf_min, sgncorr12,
                                       af0[(indp1, indp2)[daf_max_ind]],
                                       ipat, y)
                        idata += 1


    # Gate the data by spatiotemporal distance AND daf_max
    def find_bin(x, bins):
        ind = (x >= bins[:-1]) & (x < bins[1:])
        if not ind.any():
            return None
        else:
            return ind.nonzero()[0][0]

    bins_d = np.logspace(2.5, 5, 7)
    bins_af = np.array([0.1, 0.4, 0.9])
    hist2d = [[[] for b in bins_af[:-1]] for b2 in bins_d[:-1]]
    for point in data:
        ind_d = find_bin(point[0], bins_d)
        if ind_d is None:
            continue

        ind_af = find_bin(point[3], bins_af)
        if ind_af is None:
            continue

        hist2d[ind_d][ind_af].append(point)

    hist2d = np.array(hist2d, object) 
    for irow, row in enumerate(hist2d):
        for icol, cell in enumerate(row):
            hist2d[irow, icol] = np.array(cell, ndmin=2)

    if use_plot:
        fig, ax = plt.subplots()
        ax.set_xlabel('Distance [days * bp]')
        ax.set_ylabel('sgn$C_0 \cdot d\\nu_1 \cdot d\\nu_2$')
        ax.set_xscale('log')
        ax.set_xlim(0.9 * bins_d[0], 1.1 * bins_d[-1])

        # Plot different af classes separately
        for icol, col in enumerate(hist2d.T):
            ind_nonempty = [icell for icell, cell in enumerate(col) if cell.shape[1]]
            x = np.sqrt(bins_d[1:] * bins_d[:-1])[ind_nonempty]
            y = [cell[:, -1].mean() for cell in col if cell.shape[1]]
            label = 'all pats, daf  $\in$ ['+\
                    str(bins_af[icol])+', '+str(bins_af[icol+1])+'['

            ax.plot(x, y, lw=2, label=label)

        ax.legend(loc=1, fontsize=10)
        ax.grid(True)

        plt.tight_layout()
        plt.ion()
        plt.show()

    # Gate the data also by ratio days/bp for checking it's a rate process
    bins_3 = np.array([5, 20, 40])
    hist3d = [[[[]
        for b in bins_af[:-1]]
        for b2 in bins_d[:-1]]
        for b3 in bins_3[:-1]]
    for point in data:
        ind_d = find_bin(point[0], bins_d)
        if ind_d is None:
            continue

        ind_af = find_bin(point[3], bins_af)
        if ind_af is None:
            continue

        ind_3 = find_bin(1.0 * point[1] / point[2], bins_3)
        if ind_3 is None:
            continue

        hist3d[ind_3][ind_d][ind_af].append(point)

    hist3d = np.array(hist3d, object) 
    for imat, mat in enumerate(hist3d):
        for irow, row in enumerate(mat):
            for icol, cell in enumerate(row):
                hist3d[imat, irow, icol] = np.array(cell, ndmin=2)

    if use_plot:
        fig, ax = plt.subplots()
        ax.set_xlabel('Distance [days * bp]')
        ax.set_ylabel('sgn$C_0 \cdot d\\nu_1 \cdot d\\nu_2$')
        ax.set_xscale('log')
        ax.set_xlim(0.9 * bins_d[0], 1.1 * bins_d[-1])

        # Plot different af classes separately
        for imat, mat in enumerate(hist3d):
            for icol, col in enumerate(mat.T):
                ind_nonempty = [icell for icell, cell in enumerate(col) if cell.shape[1]]
                x = np.sqrt(bins_d[1:] * bins_d[:-1])[ind_nonempty]
                y = [cell[:, -1].mean() for cell in col if cell.shape[1]]
                dy = [cell[:, -1].std() for cell in col if cell.shape[1]]
                label = r'daf  $\in$ ['+\
                        str(bins_af[icol])+', '+str(bins_af[icol+1])+'[, '+\
                        '$\Delta t / \Delta pos \in$ ['+\
                        str(bins_3[imat])+', '+str(bins_3[imat+1])+'['

                ax.plot(x, y, lw=2, label=label)

        ax.legend(loc=1, fontsize=10)
        ax.grid(True)

        plt.tight_layout()
        plt.ion()
        plt.show()

