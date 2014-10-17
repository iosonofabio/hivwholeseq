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
from hivwholeseq.patients.filenames import get_initial_reference_filename
from hivwholeseq.patients.one_site_statistics import get_allele_count_trajectories



# Classes
class Propagator(object):

    def __init__(self, n_bins, use_logit=False):
        '''Prepare bins, the histogram, and other data structures'''
        if not use_logit:
            binsx = np.logspace(np.log10(0.002), np.log10(0.998), n_bins)
            binsxc = np.sqrt(binsx[1:] * binsx[:-1])
            binsy = np.concatenate([[0], binsx, [1]])
            binsyc = np.concatenate([[0], np.sqrt(binsy[2:-1] * binsy[1:-2]), [1]])
        else:
            self.trfun = trfun = lambda x: np.log10(x / (1 - x))
            self.trfuni = trfuni = lambda y: 1.0 / (1 + 10**(-y))
            binsx = trfun(np.linspace(trfuni(2e-3), trfuni(1.0 - 2e-3), n_bins))
            binsxc = trfun(0.5 * (trfuni(binsx)[1:] + trfuni(binsx)[:-1]))
            binsy = np.concatenate([[0], binsx, [1]])
            binsyc = np.concatenate([[0], trfun(0.5 * (trfuni(binsy)[2:-1] + trfuni(binsy)[1:-2])), [1]])

        binsxw = binsx[1:] - binsx[:-1]
        binsyw = binsy[1:] - binsy[:-1]

        self.binsx = binsx
        self.binsy = binsy
        self.binsxc = binsxc
        self.binsyc = binsyc
        self.binsxw = binsxw
        self.binsyw = binsyw

        self.use_logit = use_logit

        self.histogram = np.zeros((len(binsx) - 1, len(binsy) - 1), float)


    def plot(self, figaxs=None, title=''):
        '''Plot the propagator'''
        import matplotlib.pyplot as plt

        plt.ioff()

        if figaxs is None:
            fig, axs = plt.subplots(1, 2, figsize=(16, 8))
        else:
            (fig, axs) = figaxs

        # Normalize histogram
        z = self.histogram
        z = (1.0 * z.T / z.sum(axis=1)).T
        z /=  self.binsyw

        # Plot with lines
        ax = axs[0]
        for iz, zi in enumerate(z):
            xi = self.binsxc[iz]

            # Do not take the first and last final frequency bins, they include the
            # extremes (loss and fixation) and behave specially
            xf = self.binsyc[1:-1]
            y = zi[1:-1]
            xf0 = 1.2e-3
            xf1 = 1.0 - 1.2e-3

            if use_logit:
                (xf, xf0, xf1) = map(self.trfun, (xf, xf0, xf1))

            ax.plot(xf, y, lw=2, c=cm.jet(1.0 * iz / z.shape[0]),
                    label='$x_i = '+'{:1.1e}'.format(xi)+'$')
            ax.scatter(xf0, zi[0], s=80, facecolor='none', lw=2,
                       edgecolor=cm.jet(1.0 * iz / z.shape[0]))
            ax.scatter(xf1, zi[-1], s=80, facecolor='none', lw=2,
                       edgecolor=cm.jet(1.0 * iz / z.shape[0]))

        if use_logit:
            ax.set_xlim(*map(self.trfun, (1e-3, 1 - 1e-3)))
        else:
            ax.set_xscale('log')
            ax.set_xlim(1e-3, 1.5)
        
        ax.grid(True)
        ax.set_ylabel('P(x1 | x0)')
        ax.set_xlabel('Final frequency')
        ax.set_yscale('log')
        ax.legend(loc=3, fontsize=10, title='Initial frequency:', ncol=2)

        # Plot with heatmap
        ax = axs[1]

        # Do not take the first and last final frequency bins, they include the
        # extremes (loss and fixation) and behave specially
        z1 = np.log10(z[:, 1:-1])

        im = ax.imshow(z1.T, interpolation='nearest')
        ax.set_xlabel('Initial freq')
        ax.set_ylabel('Final freq')
        ax.set_xticks(np.arange(len(self.binsx)) - 0.5)
        ax.set_xticklabels(map('{:1.2e}'.format, self.binsx),
                           rotation=45, fontsize=10)
        ax.set_yticks(np.arange(len(self.binsy) - 2) - 0.5)
        ax.set_yticklabels(map('{:1.2e}'.format, self.binsy[1:-1]),
                           rotation=45, fontsize=10)

        # Reverse the y axis (by default image y coordinates are top to bottom)
        ax.set_ylim(*(ax.get_ylim()[::-1]))

        cb = plt.colorbar(im)
        cb.set_label('log10 P(x1 | x0)', labelpad=30, rotation=270, fontsize=12)

        if title:
            fig.suptitle(title, fontsize=16)
        plt.tight_layout(rect=(0, 0, 1, 0.94))



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

    # Prepare output structures
    n_bins = 14
    pp = Propagator(n_bins, use_logit=use_logit)

    binsd = np.array([-0.5, -0.4, -0.3, -0.2, -0.1, -0.05, -0.02, -0.01, -0.003,
                      0.003, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5])
    binsdc = 0.5 * (binsd[1:] + binsd[:-1])
    binsdw = binsd[1:] - binsd[:-1]

    histr = pp.histogram
    histd = np.zeros((len(pp.binsx) - 1, len(binsd) - 1), float)

    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        samplenames = patient.samples.index

        # If the script is called with no fragment, iterate over all
        if not fragments:
            fragments = ['F'+str(i) for i in xrange(1, 7)]
        if VERBOSE >= 2:
            print 'fragments', fragments
    
        # Iterate over samples and fragments
        for fragment in fragments:
            if VERBOSE >= 1:
                print pname, fragment
    
            aft, ind = patient.get_allele_frequency_trajectories(fragment,
                                                                 cov_min=depth_min)

            n_templates = np.array(patient.n_templates[ind])
            indd = n_templates >= depth_min
            aft = aft[indd]
            ind = ind[indd]
            n_templates = n_templates[indd]

            ts = patient.times[ind]
    
            # Collect counts
            for i in xrange(aft.shape[0] - 1):
                for j in xrange(i + 1, aft.shape[0]):
                    dij = j - i
                    if ((ts[j] - ts[i]) > dt[1]) or ((ts[j] - ts[i]) < dt[0]):
                        continue

                    histr += np.histogram2d(aft[i].ravel(),
                                            aft[j].ravel(),
                                            bins=[pp.binsx, pp.binsy])[0]

                    histd += np.histogram2d(aft[i].ravel(),
                                            aft[j].ravel() - aft[i].ravel(),
                                            bins=[pp.binsx, binsd])[0]

    if plot:
        pp.plot(title='Propagator for allele frequencies\n$\Delta t = '+str(dt)+'$ days, '+str(fragments))
        plt.ion()
        plt.show()

        sys.exit()

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
        ax.legend(loc=1, fontsize=10, title='Initial\nfrequency:', ncol=2)

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
        ax.legend(loc=3, fontsize=10, title='Initial\nfrequency:', ncol=2)

        ax.set_title('Propagator (normalized )for allele frequencies\n'+\
                     '$\Delta t = '+str(dt)+'$, '+str(fragments),
                     fontsize=16)
        plt.tight_layout()


        plt.ion()
        plt.show()
