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

from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.argparse_utils import RoiAction



# Functions
def plot_logo(af,
              alphabet=('A', 'C', 'G', 'T'),
              colors=('green', 'blue', 'orange', 'red'),
              ax=None, VERBOSE=0,
              xoffset=0, yoffset=0):
    '''Print a logo plot (entropy style)'''
    from matplotlib.patches import Rectangle as R

    S = (-1.0 * af * np.log2(af + 1e-10)).sum(axis=0)
    info = 2.0 - S
    heights = af * info
    ax.set_xlim(-1 + xoffset, heights.shape[1] + 1 + xoffset)
    ax.set_ylim(-0.05 + yoffset, 2.05 + yoffset)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    for pos, hs in enumerate(heights.T):
        if VERBOSE >= 3:
            print pos, hs

        hs[hs < 1e-2] = 0
        hs[hs > 2 - 1e-2] = 2
        if (hs == 0).all():
            continue

        #import ipdb; ipdb.set_trace()

        inds = np.argsort(hs)[::-1]
        hs = hs[inds]
        hscum = np.concatenate([[0], hs.cumsum()])
        for ii, i in enumerate(inds):
            h = hs[ii]
            if h > 0:
                y0 = hscum[ii]
                y1 = hscum[ii + 1]
                rect = R((pos + xoffset, y0 + yoffset), width=1, height=(y1 - y0),
                         facecolor=colors[i], edgecolor='k', lw=0.5,
                         alpha=0.6)
                ax.add_patch(rect)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get interpatient diversity',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--roi', action=RoiAction, required=True,
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot the allele frequency trajectories')

    args = parser.parse_args()
    pnames = args.patients
    roi = args.roi
    VERBOSE = args.verbose
    use_plot = args.plot

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]
    
    refname = 'HXB2'
    (fragment, roi_start, roi_end) = roi

    # Collect allele counts
    alleles_inter = np.zeros((4, 10000), int)
    af_intra = np.ma.masked_all((len(patients), 4, 10000), float)
    frag_start = 100000000
    frag_end = 0
    for ip, (pname, patient) in enumerate(patients.iterrows()):
        patient = Patient(patient)

        mapco = patient.get_map_coordinates_reference(fragment, refname=refname)
        frag_start = min(frag_start, mapco[:, 0].min())
        frag_end = max(frag_end, mapco[:, 0].max())

        act, ind = patient.get_allele_count_trajectories(fragment)
        for (pos_ref, pos_pat) in mapco:
            # Interpatient
            abu = Counter(act[:, :4, pos_pat].argmax(axis=1))
            for i, count in abu.iteritems():
                alleles_inter[i, pos_ref] += count

            # Intrapatient
            abu = act[:, :4, pos_pat]
            if (abu.sum(axis=1) > 10).all():
                abu = (1.0 * abu.T / abu.sum(axis=1)).T
                af_intra[ip, :, pos_ref] = abu.mean(axis=0)

    # Interpatient
    alleles_inter = alleles_inter[:, frag_start + roi_start: frag_start + roi_end]
    af_inter = np.ma.array(1.0 * alleles_inter / (1e-3 + alleles_inter.sum(axis=0)),
                           mask=np.tile((alleles_inter == 0).all(axis=0), [4, 1]),
                           dtype=float)

    if use_plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(23, 4))
        alphabet=('A', 'C', 'G', 'T')
        colors=('green', 'blue', 'orange', 'red')
        plot_logo(af_inter, ax=ax, alphabet=alphabet, colors=colors, xoffset=roi_start)
        ax.set_title('Sequence logo across patients, '+fragment)
        ax.set_xlabel('Position [bp]')
        ax.set_ylabel('Shannon information [bits]')
        plt.tight_layout()

    # Intrapatient
    af_intra = af_intra[:, :, frag_start + roi_start: frag_start + roi_end]

    if use_plot:
        import matplotlib.pyplot as plt
        alphabet=('A', 'C', 'G', 'T')
        colors=('green', 'blue', 'orange', 'red')
        for ip, af_intra_pat in enumerate(af_intra):
            fig, ax = plt.subplots(figsize=(23, 4))
            ax.set_title('Sequence logo, patient '+patients.iloc[ip].name+', '+fragment)
            plot_logo(af_intra_pat, ax=ax, alphabet=alphabet, colors=colors,
                      xoffset=roi_start)
            ax.set_ylabel('Shannon information [bits]')
            ax.set_xlabel('Position [bp]')
            plt.tight_layout()

    if use_plot:
        plt.ion()
        plt.show()
