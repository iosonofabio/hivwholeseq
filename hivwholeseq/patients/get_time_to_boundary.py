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

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.patients import load_patient



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get time to boundary (fix/loss)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 genomewide)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot the time distributions')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose
    plot = args.plot
    cov_min = 200
    depth_min = 100

    patient = load_patient(pname)
    patient.discard_nonsequenced_samples()

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    af0 = [0.15, 0.85]
    af_bd = [0.01, 0.99]
    t_bds = []
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

        t_bd = []
        for pos in xrange(aft.shape[2]):
            for ia, a in enumerate(alpha):
                aft_pos = aft[:, ia, pos]
                
                # Keep only polymorphic
                ipos0 = (aft_pos > af0[0]) & (aft_pos < af0[1])
                if not ipos0.any():
                    continue

                it0 = ipos0.nonzero()[0][0]

                # Keep only if they fix/extinct within temporal window
                iposbd = (aft_pos[it0:] < af_bd[0]) | (aft_pos[it0:] > af_bd[1])
                if not iposbd.any():
                    continue

                itbd = it0 + iposbd.nonzero()[0][0]
                t_bd.append(times[itbd] - times[it0])


        t_bds.append(t_bd)

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
        ax.set_title('Time to boundary from '+'-'.join(map(str, af0))\
                     +', patient '+patient.name)

        plt.ion()
        plt.show()

