# vim: fdm=marker
'''
author:     Fabio Zanini
date:       02/07/14
content:    Explore parallel evolution in different patients by studying shared
            allele frequency trajectories.
'''
# Modules
import os
import argparse
from operator import itemgetter
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO, AlignIO




# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get shared allele trajectories',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')

    args = parser.parse_args()
    fragments = args.fragments
    VERBOSE = args.verbose

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    plt.ioff()

    for fragment in fragments:
        if VERBOSE >= 1:
            print fragment

        # Low-coverage regions are bytecoded by -1
        from hivwholeseq.patients.filenames import root_patient_folder
        afts = np.load(root_patient_folder+'all/aft_shared_'+fragment+'.npy')
        afts = afts.swapaxes(0, 1).swapaxes(1, 2)
        n_patients = afts.shape[2]

        # Look for parallel positive selection
        hist_sweep = np.zeros(n_patients + 1, int)
        for pos, aft_pos in enumerate(afts):
            if VERBOSE >= 3:
                print pos
            for aft_nuc in aft_pos:
                n = 0
                for aft_pat in aft_nuc:
                    if (-0.1 < aft_pat[0] < 0.5) and (aft_pat[1:] > 0.5).any():
                        n += 1
                hist_sweep[n] += 1

        print '{:25s}'.format('Sweep sharing:')+' '.join(map('{:>7d}'.format, hist_sweep))

        # Look for parallel conservation
        threshold = 0.005
        hist_cons = np.zeros(n_patients + 1, int)
        for pos, aft_pos in enumerate(afts):
            if VERBOSE >= 3:
                print pos
            for aft_nuc in aft_pos:
                n = 0
                for aft_pat in aft_nuc:
                    # FIXME: deal better with missing data (low coverage)
                    times_cov = aft_pat >= -0.1
                    if times_cov.sum() < 2:
                        continue
                    aft_pat = aft_pat[times_cov]
                    if ((aft_pat > threshold) & (aft_pat < 1 - threshold)).any():
                        n += 1
                hist_cons[n_patients - n] += 1

        print '{:25s}'.format('Cons sharing ('+str(threshold * 100)+'%):')+' '.join(map('{:>7d}'.format, hist_cons))

        # Build measure of conservation based on that
        n_cov = np.zeros(afts.shape[0], int)
        n_poly = np.zeros((afts.shape[0], afts.shape[1]), int)
        for pos, aft_pos in enumerate(afts):
            if VERBOSE >= 3:
                print pos
            for j, aft_nuc in enumerate(aft_pos):
                n = 0
                nc = 0
                for aft_pat in aft_nuc:
                    # FIXME: deal better with missing data (low coverage)
                    times_cov = aft_pat >= -0.1
                    if times_cov.sum() < 2:
                        continue
                    nc += 1
                    aft_pat = aft_pat[times_cov]
                    if ((aft_pat > threshold) & (aft_pat < 1 - threshold)).any():
                        n += 1
                n_cov[pos] = nc
                n_poly[pos, j] = n

        n_poly_max = n_poly.max(axis=1)
        n_poly_max = np.ma.array(n_poly_max, mask=(n_cov < 3))
        fig, ax = plt.subplots(figsize=(16, 5))
        ax.plot(np.arange(len(n_poly_max)), n_poly_max, lw=1.5)
        ax.set_xlabel('Position [bp]')
        ax.set_ylabel('# patients polymorphic')
        ax.set_title(fragment)
        ax.set_xlim(-1, len(n_poly_max) + 1)

    plt.tight_layout()
    plt.ion()
    plt.show()

