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

from hivwholeseq.miseq import alphal
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.argparse_utils import RoiAction



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
    parser.add_argument('--plot', default='both',
                        choices=['intra', 'inter', 'cross', 'pats', 'intercross', 'all'],
                        help='Plot the allele frequency trajectories')
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
    
    refname = 'HXB2'
    (fragment, roi_start, roi_end) = roi

    # Make patient classes (FIXME: do it better?)
    patients = [Patient(p) for pn, p in patients.iterrows()]

    # Get coordinates
    for patient in patients:
        patient.mapco = patient.get_map_coordinates_reference(fragment, refname=refname)

    # Get fragment coords
    frag_start = min([patient.mapco[:, 0].min() for patient in patients])
    frag_end = max([patient.mapco[:, 0].max() for patient in patients])

    # Collect allele counts
    alleles_inter = np.zeros((4, roi_end - roi_start), int)
    af_intra = np.ma.masked_all((len(patients), 4, roi_end - roi_start), float)
    for ip, patient in enumerate(patients):

        mapco = patient.mapco
        mapco = mapco[(mapco[:, 0] >= frag_start + roi_start) & \
                      (mapco[:, 0] < frag_start + roi_end)]

        act, ind = patient.get_allele_count_trajectories(fragment)
        for (pos_ref, pos_pat) in mapco:
            # Interpatient
            abu = Counter(act[:, :4, pos_pat].argmax(axis=1))
            for i, count in abu.iteritems():
                alleles_inter[i, pos_ref - frag_start - roi_start] += count
            # TODO: decide on gaps!

            # Intrapatient
            abu = act[:, :4, pos_pat]
            if (abu.sum(axis=1) > 10).all():
                abu = (1.0 * abu.T / abu.sum(axis=1)).T
                af_intra[ip, :, pos_ref - frag_start - roi_start] = abu.mean(axis=0)

    # Interpatient
    af_inter = np.ma.array(1.0 * alleles_inter / (1e-3 + alleles_inter.sum(axis=0)),
                           mask=np.tile((alleles_inter == 0).all(axis=0), [4, 1]),
                           dtype=float)

    if use_plot in ['inter', 'pats', 'all', 'intercross']:
        fig, ax = plt.subplots(figsize=(23, 4))
        alphabet=('A', 'C', 'G', 'T')
        colors=('green', 'blue', 'orange', 'red')
        plot_logo(af_inter, ax=ax, alphabet=alphabet, colors=colors,
                  xoffset=frag_start + roi_start, picker=use_interactive)
        ax.set_title('Sequence logo across patients, '+fragment)
        ax.set_xlabel('Position [bp]')
        ax.set_ylabel('Shannon information [bits]')
        plt.tight_layout()

        if use_interactive:
            fig.canvas.mpl_connect('pick_event', onpick_callback)


    # Intrapatient
    if use_plot in ['intra', 'pats', 'all']:
        alphabet=('A', 'C', 'G', 'T')
        colors=('green', 'blue', 'orange', 'red')
        for ip, af_intra_pat in enumerate(af_intra):
            fig, ax = plt.subplots(figsize=(23, 4))
            ax.set_title('Sequence logo, patient '+patients[ip].name+', '+fragment)
            plot_logo(af_intra_pat, ax=ax, alphabet=alphabet, colors=colors,
                      xoffset=frag_start + roi_start, picker=use_interactive)
            ax.set_ylabel('Shannon information [bits]')
            ax.set_xlabel('Position [bp]')
            plt.tight_layout()

            if use_interactive:
                fig.canvas.mpl_connect('pick_event', onpick_callback)

    if use_plot:
        plt.ion()
        plt.show()


    # Cross-sectional
    import glob
    from Bio import AlignIO
    csdata_folder = '/ebio/ag-neher/home/fzanini/phd/general/data/typeM/'
    from hivwholeseq.reference import load_custom_reference
    refseq = load_custom_reference(refname, 'gb')
    # Check all genes
    for fea in refseq.features:
        feaname = fea.id

        # Check for the MSA file
        fns = glob.glob(csdata_folder+'*'+feaname+'*'+'.fasta')
        if not len(fns):
            continue
        fn = fns[0]

        # Check for overlap with our roi (assume only one is good)
        # FIXME: deal with exons
        if feaname in ('tat', 'rev'):
            continue

        fea_start = fea.location.nofuzzy_start
        fea_end = fea.location.nofuzzy_end
        # Full coverage
        if (fea_start <= frag_start + roi_start) and (fea_end >= frag_start + roi_end):
            ali = AlignIO.read(fn, 'fasta')
            break

        # TODO: implement partial coverage

    else:
        sys.exit()

    # Delete alignment columns that are absent from the chosen reference 
    # for now assume the reference is part of the ali, it's much faster that way
    for seq in ali:
        if refname in seq.id:
            break
    else:
        raise ValueError('Reference not found in the alignment')

    mapcoref = []
    pos_ref = fea_start
    for pos, nuc in enumerate(seq.seq):
        if nuc != '-':
            mapcoref.append((pos_ref, pos))
            pos_ref += 1
    mapcoref = np.array(mapcoref, int)
    mapcoref = mapcoref[(mapcoref[:, 0] >= frag_start + roi_start) & \
                        (mapcoref[:, 0] < frag_start + roi_end)]
 
    af_cross = np.ma.masked_all((5, roi_end - roi_start), float)
    # Gaps?
    alpha4 = alphal[:4]
    for (pos_ref, pos_ali) in mapcoref:
        alleles = Counter(ali[:, pos_ali])
        for i, a in enumerate(alpha4):
            af_cross[i, pos_ref - frag_start - roi_start] = alleles[a]
            
        af_cross[:, pos_ref - frag_start - roi_start] /= \
                af_cross[:, pos_ref - frag_start - roi_start].sum()

    if use_plot in ('all', 'cross', 'intercross'):
        # Gaps?
        alphabet=tuple(alpha4)
        colors=('green', 'blue', 'orange', 'red')
        fig, ax = plt.subplots(figsize=(23, 4))
        ax.set_title('Sequence logo, cross-sectional (type M)')
        plot_logo(af_cross, ax=ax, alphabet=alphabet, colors=colors,
                  xoffset=frag_start + roi_start, picker=use_interactive)
        ax.set_ylabel('Shannon information [bits]')
        ax.set_xlabel('Position [bp]')
        plt.tight_layout()

        if use_interactive:
            fig.canvas.mpl_connect('pick_event', onpick_callback)

    if use_plot:
        plt.ion()
        plt.show()

