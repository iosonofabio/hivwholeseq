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

from hivwholeseq.patients.patients import load_patients, filter_patients_n_times, Patient
from hivwholeseq.patients.filenames import get_initial_consensus_filename
from hivwholeseq.patients.one_site_statistics import get_allele_count_trajectories



# Script
if __name__ == '__main__': 

    # Parse input args
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
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the propagator')
    parser.add_argument('--deltat', type=int, default=1,
                        help='Number of time points between final and initial')
    parser.add_argument('--logit', action='store_true',
                        help='use logit scale (log(x/(1-x)) in the plots')


    args = parser.parse_args()
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    save_to_file = args.save
    plot = args.plot
    submit = args.submit
    dt = args.deltat
    use_logit = args.logit

    # Restrict to patients with at least three time points
    patients = load_patients()
    ind = filter_patients_n_times(patients, n_times=3)
    ind &= patients['Initial reference'].notnull()
    ind &= -patients.index.isin(['15034'])  # FIXME: Superinfection? Mislabelling?
    if pnames is not None:
        ind &= patients.index.isin(pnames)
    patients = patients.loc[ind]

    # Prepare output structures
    n_bins = 14
    if not use_logit:
        binsx = np.logspace(np.log10(0.002), np.log10(0.998), n_bins)
        binsxc = np.sqrt(binsx[1:] * binsx[:-1])
        binsy = np.concatenate([[0], binsx, [1]])
        binsyc = np.concatenate([[0], np.sqrt(binsy[2:-1] * binsy[1:-2]), [1]])
    else:
        trfun = lambda x: np.log10(x / (1 - x))
        trfuni = lambda y: 1.0 / (1 + 10**(-y))
        binsx = trfun(np.linspace(trfuni(2e-3), trfuni(1.0 - 2e-3), n_bins))
        binsxc = trfun(0.5 * (trfuni(binsx)[1:] + trfuni(binsx)[:-1]))
        binsy = np.concatenate([[0], binsx, [1]])
        binsyc = np.concatenate([[0], trfun(0.5 * (trfuni(binsy)[2:-1] + trfuni(binsy)[1:-2])), [1]])

    binsxw = binsx[1:] - binsx[:-1]
    binsyw = binsy[1:] - binsy[:-1]

    binsd = np.array([-0.5, -0.4, -0.3, -0.2, -0.1, -0.05, -0.02, -0.01, -0.003,
                      0.003, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5])
    binsdc = 0.5 * (binsd[1:] + binsd[:-1])
    binsdw = binsd[1:] - binsd[:-1]

    histr = np.zeros((len(binsx) - 1, len(binsy) - 1), float)
    histd = np.zeros((len(binsx) - 1, len(binsd) - 1), float)

    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()
        samplenames = patient.samples.index

        # If the script is called with no fragment, iterate over all
        if not fragments:
            fragments = ['F'+str(i) for i in xrange(1, 7)]
        if VERBOSE >= 2:
            print 'fragments', fragments
    
        # Iterate over samples and fragments
        for fragment in fragments:
    
            # Submit to the cluster self if requested (--save is assumed)
            if submit:
                fork_self(pname, fragment,
                          VERBOSE=VERBOSE)
                continue
    
            if VERBOSE >= 1:
                print pname, fragment
    
            # Collect allele counts from patient samples, and return only positive hits
            # sns contains sample names and PCR types
            # Keep PCR2 only if PCR1 is absent
            (sns, act) = get_allele_count_trajectories(pname, samplenames, fragment,
                                                       use_PCR1=1, VERBOSE=VERBOSE)

            if len(sns) < 1 + dt:
                continue

            ind = [i for i, (_, sample) in enumerate(patient.samples.iterrows())
                   if sample.name in map(itemgetter(0), sns)]
            samples = patient.samples.iloc[ind]
            #times = (samples.date - patient.transmission_date) / np.timedelta64(1, 'D')
            ntemplates = samples['n templates']
            depthmax = np.maximum(1.0 / np.array(ntemplates), 1e-3)

            aft = (1.0 * act.swapaxes(0, 1) / act.sum(axis=1)).swapaxes(0, 1)
            aft[(aft.swapaxes(0, 2) < depthmax).swapaxes(0, 2)] = 0
            aft[(aft.swapaxes(0, 2) > 1 - depthmax).swapaxes(0, 2)] = 1
            aft = (aft.swapaxes(0, 1) / aft.sum(axis=1)).swapaxes(0, 1)
    
            # Collect counts
            for i in xrange(aft.shape[0] - dt):
                histr += np.histogram2d(aft[i].ravel(),
                                        aft[i + dt].ravel(),
                                        bins=[binsx, binsy])[0]

                histd += np.histogram2d(aft[i].ravel(),
                                        aft[i + dt].ravel() - aft[i].ravel(),
                                        bins=[binsx, binsd])[0]

    if plot:
        plt.ioff()

        fig, axs = plt.subplots(1, 2, figsize=(16, 8))
        z = histr
        z = (1.0 * z.T / z.sum(axis=1)).T
        z /=  binsyw

        ax = axs[0]
        for iz, zi in enumerate(z):
            xi = binsxc[iz]
            xf = binsyc[1:-1]
            y = zi[1:-1]
            xf0 = 1.2e-3
            xf1 = 1.0 - 1.2e-3

            if use_logit:
                (xf, xf0, xf1) = map(trfun, (xf, xf0, xf1))

            ax.plot(xf, y, lw=2, c=cm.jet(1.0 * iz / z.shape[0]),
                    label='$x_i = '+'{:1.1e}'.format(xi)+'$')
            ax.scatter(xf0, zi[0], s=80, facecolor='none', lw=2,
                       edgecolor=cm.jet(1.0 * iz / z.shape[0]))
            ax.scatter(xf1, zi[-1], s=80, facecolor='none', lw=2,
                       edgecolor=cm.jet(1.0 * iz / z.shape[0]))

        if use_logit:
            ax.set_xlim(*map(trfun, (1e-3, 1 - 1e-3)))
        else:
            ax.set_xscale('log')
            ax.set_xlim(1e-3, 1.5)

        ax.set_ylabel('P(x1 | x0)')
        ax.set_xlabel('Final frequency')
        ax.set_yscale('log')
        ax.legend(loc=3, fontsize=10)

        ax = axs[1]
        z1 = np.log10(z[:, 1:-1])
        im = ax.imshow(z1.T, interpolation='nearest')
        ax.set_xlabel('Initial freq')
        ax.set_ylabel('Final freq')
        ax.set_xticks(np.arange(len(binsx)) - 0.5)
        ax.set_xticklabels(map('{:1.1e}'.format, binsx),
                           rotation=45, fontsize=10)
        ax.set_yticks(np.arange(len(binsy) - 2) - 0.5)
        ax.set_yticklabels(map('{:1.2e}'.format, binsy[1:-1]),
                           rotation=45, fontsize=10)
        ax.set_ylim(*(ax.get_ylim()[::-1]))
        cb = plt.colorbar(im)
        cb.set_label('log10 P(x1 | x0)', labelpad=30, rotation=270, fontsize=12)

        fig.suptitle('Propagator for allele frequencies\n'+\
                     '$\Delta t = '+str(dt)+'$, '+str(fragments),
                     fontsize=16)
        plt.tight_layout(rect=(0, 0, 1, 0.94))

        # Plot difference
        fig, axs = plt.subplots(1, 2, figsize=(16, 8))
        z = histd
        z = 1.0 * (z.T / z.sum(axis=1)).T
        z /= binsdw

        ax = axs[0]
        for iz, zi in enumerate(z):
            xi = binsxc[iz]
            dx = binsdc[1:-1]
            y = zi[1:-1]

            ax.plot(dx, y, lw=2, c=cm.jet(1.0 * iz / z.shape[0]),
                    label='$x_i = '+'{:1.1e}'.format(xi)+'$')

        ax.set_ylabel('P(x1 - x0 | x0)')
        ax.set_xlabel('x1 - x0')
        ax.set_yscale('log')
        ax.legend(loc=1, fontsize=10)

        ax = axs[1]
        z1 = np.log10(z)
        z1 = np.maximum(z1, z1[z1 != (-np.inf)].min() - 1)
        im = ax.imshow(z1.T, interpolation='nearest', aspect='auto')
        ax.set_xlabel('Initial freq')
        ax.set_ylabel('Freq diff')
        ax.set_xticks(np.arange(len(binsx)) - 0.5)
        ax.set_xticklabels(map('{:1.1e}'.format, binsx),
                           rotation=45, fontsize=10)
        ax.set_yticks(np.arange(len(binsd)) - 0.5)
        ax.set_yticklabels(map('{:1.2e}'.format, binsd),
                           rotation=45, fontsize=10)
        ax.set_ylim(*(ax.get_ylim()[::-1]))
        cb = plt.colorbar(im)
        cb.set_label('log10 P(x1 - x0 | x0)', labelpad=30, rotation=270, fontsize=12)

        fig.suptitle('Propagator for allele frequencies\n'+\
                     '$\Delta t = '+str(dt)+'$, '+str(fragments),
                     fontsize=16)
        plt.tight_layout(rect=(0, 0, 1, 0.94))

        # Plot normalized (Kosheleva et Desai 2013)
        # rho(x_k | x_k-1) = x_k-1 (1 - x_k-1) / (q * dx^2) -- q = 8 in HIV
        fig, ax = plt.subplots(1)
        z = histr
        z = (1.0 * z.T / z.sum(axis=1)).T
        z /=  binsyw

        for iz, zi in enumerate(z):
            xi = binsxc[iz]
            xf = binsyc[1:-1]
            y = zi[1:-1] / (xi * (1 - xi)) * ((xf - xi)**2)

            if use_logit:
                (xf, xf0, xf1) = map(trfun, (xf, xf0, xf1))

            ax.plot(xf, y, lw=2, c=cm.jet(1.0 * iz / z.shape[0]),
                    label='$x_i = '+'{:1.1e}'.format(xi)+'$')

        if use_logit:
            ax.set_xlim(*map(trfun, (2e-3, 1 - 2e-3)))
        else:
            ax.set_xscale('log')
            ax.set_xlim(1e-3, 1.5)

        ax.set_ylabel('P(x1 | x0) / (x0 * (1 - x0)) * dx^2')
        ax.set_xlabel('Final frequency')
        ax.set_yscale('log')
        ax.legend(loc=3, fontsize=10)

        ax.set_title('Propagator (normalized )for allele frequencies\n'+\
                     '$\Delta t = '+str(dt)+'$, '+str(fragments),
                     fontsize=16)
        plt.tight_layout()



        plt.ion()
        plt.show()
