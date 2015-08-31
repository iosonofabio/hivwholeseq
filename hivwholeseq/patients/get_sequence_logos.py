# vim: fdm=marker
'''
author:     Fabio Zanini
date:       15/10/14
content:    Get interpatient diversity along the genome.
'''
# Modules
import os, sys
import argparse
from itertools import izip
from collections import defaultdict, Counter
import numpy as np
import matplotlib.pyplot as plt

from hivwholeseq.utils.miseq import alphal
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.utils.argparse import RoiAction



# Functions
def plot_logo(af,
              alphabet=('A', 'C', 'G', 'T'),
              colors=('green', 'blue', 'orange', 'red'),
              ax=None, VERBOSE=0,
              xoffset=0, yoffset=0,
              reverse=True,
              **kwargs):
    '''Print a logo plot (Shannon info).
    
    Parameters:
      **kwargs: passed to matplotlib.patches.Rectangle
    '''
    import numpy as np
    from matplotlib.patches import Rectangle as R

    S = (-1.0 * af * np.log2(af + 1e-10)).sum(axis=0)
    Smax = np.log2(len(alphabet))
    info = Smax - S
    heights = af * info
    ax.set_xlim(-1.5 + xoffset, heights.shape[1] + 0.5 + xoffset)
    ax.set_ylim(-0.05 + yoffset, Smax + 0.05 + yoffset)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    for pos, hs in enumerate(heights.T):
        if VERBOSE >= 3:
            print pos, hs

        hs[hs < 1e-2] = 0
        hs[hs > 2 - 1e-2] = 2
        if (hs == 0).all():
            continue

        inds = np.argsort(hs)
        if reverse:
            inds = inds[::-1]
        hs = hs[inds]
        hscum = np.concatenate([[0], hs.cumsum()])
        for ii, i in enumerate(inds):
            h = hs[ii]
            if h > 0:
                y0 = hscum[ii]
                y1 = hscum[ii + 1]
                rect = R((-0.5 + pos + xoffset, y0 + yoffset), width=1, height=(y1 - y0),
                         facecolor=colors[i], edgecolor='k', lw=0.5,
                         label=alphabet[i],
                         alpha=0.6, **kwargs)
                ax.add_patch(rect)


def onpick_callback(event):
    '''Print the base corresponding to a rectangle'''
    rect = event.artist
    nuc = rect.get_label()
    print nuc



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Get interpatient diversity',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--roi', action=RoiAction, required=True,
                        help='Region of Interest (e.g. F1 100 200)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot the logos')
    parser.add_argument('--interactive', action='store_true',
                        help='Add mouse events to the plot')

    args = parser.parse_args()
    pnames = args.patients
    roi = args.roi
    VERBOSE = args.verbose
    use_plot = args.plot
    use_interactive = args.interactive

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]
    

    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        (fragment, start, end) = patient.get_fragmented_roi(roi, VERBOSE=VERBOSE)
        aft, ind = patient.get_allele_frequency_trajectories(fragment)
        aft = aft[:, :, start: end]

        # TODO: also calculate the logos

        ## Get only some time points
        #i = np.arange(len(ind))[::len(ind) // 2]
        #aft = aft[i]
        #ind = ind[i]

        times = patient.times[ind]

        if use_plot:
            fig, axs = plt.subplots(aft.shape[0], 1, figsize=(14, 3 * aft.shape[0]))

            for i, (ax, af) in enumerate(izip(axs, aft)):
                plot_logo(af[:4], ax=ax, VERBOSE=VERBOSE,
                          xoffset=roi[1])
                ax.set_ylabel(str(times[i])+' days')
                if i < len(axs) - 1:
                    ax.set_xticks([])
                else:
                    ax.set_xlabel('Position in '+roi[0]+' [bp]')

            plt.tight_layout(h_pad=0.001)

    if use_plot:
        plt.ion()
        plt.show()
