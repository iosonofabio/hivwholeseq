# vim: fdm=marker
'''
author:     Fabio Zanini
date:       12/11/14
content:    Measure the time for polymorphisms to reach the boundary (fixation
            or loss), we might use it to estimate coalescence time.
'''
# Modules
import os
import argparse
from itertools import izip
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.utils.miseq import alpha
from hivwholeseq.patients.patients import load_patients, Patient



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get time to boundary (fix/loss)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 genomewide)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot the time distributions')

    args = parser.parse_args()
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    plot = args.plot
    cov_min = 200
    depth_min = 100

    af0 = [0.15, 0.85]
    af_bd = [0.05, 0.95]

    data = {}

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    for pname, patient in patients.iterrows():
        print pname
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()

        t_bds = []
        t_loss = []
        t_fixs = []
        n_staypolys = []
        for fragment in fragments:
            if VERBOSE >= 1:
                print fragment

            # Collect allele counts from patient samples, and return only positive hits
            # sns contains sample names and PCR types
            (aft, ind) = patient.get_allele_frequency_trajectories(fragment,
                                                               cov_min=cov_min,
                                                               depth_min=depth_min,
                                                               VERBOSE=VERBOSE)
            times = patient.times[ind]
            ntemplates = patient.n_templates[ind]

            n_staypoly = 0
            t_bd = []
            t_fix = []
            t_los = []
            for pos in xrange(aft.shape[2]):
                for ia, a in enumerate(alpha):
                    aft_pos = aft[:, ia, pos]
                    
                    # Keep only polymorphic
                    ipos0 = (aft_pos > af0[0]) & (aft_pos < af0[1])
                    if not ipos0.any():
                        continue

                    it0 = ipos0.nonzero()[0][0]

                    # If they do not fix/extinct within temporal window, assign
                    # long time
                    iposbd = (aft_pos[it0:] < af_bd[0]) | (aft_pos[it0:] > af_bd[1])
                    if not iposbd.any():
                        t_bd.append(10000)
                        n_staypoly += 1
                        continue

                    itbd = it0 + iposbd.nonzero()[0][0]
                    t_bd.append(times[itbd] - times[it0])

                    iposfix = (aft_pos[it0:] > af_bd[1])
                    iposlos = (aft_pos[it0:] < af_bd[0])
                    if iposfix.any() and (not iposlos.any()):
                        t_fix.append(times[it0 + iposfix.nonzero()[0][0]] - times[it0])
                    elif (not iposfix.any()) and iposlos.any():
                        t_los.append(times[it0 + iposlos.nonzero()[0][0]] - times[it0])
                    else:
                        itfix = it0 + iposfix.nonzero()[0][0]
                        itlos = it0 + iposlos.nonzero()[0][0]
                        # It cannot be both fixed and lost at the same time
                        if itfix < itlos:
                            t_fix.append(times[itfix] - times[it0])
                        else:
                            t_los.append(times[itlos] - times[it0])


            t_bds.append(t_bd)
            t_fixs.append(t_fix)
            t_loss.append(t_los)
            n_staypolys.append(n_staypoly)

        data[pname] = {'t_boundary': t_bds, 't_fix': t_fixs, 't_los': t_loss,
                       'staypoly': n_staypolys}

        if plot:
            fig, ax = plt.subplots()
            for ifr, (fragment, t_bd) in enumerate(izip(fragments, t_bds)):
                x = np.sort(t_bd)
                y = np.linspace(0, 1, len(x))[::-1]
                ax.plot(x, y, label=fragment, lw=2,
                        color=cm.jet(1.0 * ifr / len(fragments)))

            ax.set_xlabel('Time to boundary [days]')
            ax.set_ylabel('Cumulative distribution')
            ax.grid(True)
            ax.legend(loc=1, fontsize=12, title='Fragments')
            ax.set_title('Time to boundary from '+'-'.join(map(str, af0))+\
                         ', patient '+patient.name)
            ax.set_xlim(xmax=x[x<10000].max()+1)
            ax.set_ylim(-0.01, 1.01)

            fig, ax = plt.subplots()
            for ifr, fragment in enumerate(fragments):
                xfix = np.sort(t_fixs[ifr])
                xlos = np.sort(t_loss[ifr])
                n_staypoly = n_staypolys[ifr]
                n_tot = n_staypoly + len(xfix) + len(xlos)
                yfix = 1 - np.linspace(0, 1.0 * len(xfix) / n_tot, len(xfix))
                ylos = np.linspace(0, 1.0 * len(xlos) / n_tot, len(xlos))
                ax.plot(xfix, yfix, lw=2,
                        label=fragment,
                        color=cm.jet(1.0 * ifr / len(fragments)))
                ax.plot(xlos, ylos, lw=2,
                        color=cm.jet(1.0 * ifr / len(fragments)))

            ax.set_xlabel('Time to boundary [days]')
            ax.set_ylabel('Cumulative distribution')
            ax.grid(True)
            ax.legend(loc=1, fontsize=12, title='Fragments')
            ax.set_title('Time to boundary from '+'-'.join(map(str, af0))+\
                         ', patient '+patient.name)
            ax.set_xlim(xmax=max(xfix.max(), xlos.max())+1)
            ax.set_ylim(-0.01, 1.01)


            plt.ion()
            plt.show()

