# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/10/14
content:    Study whether HIV in a patient tends to explore the mutational space
            in a way that comes closer to a subtype average (entropy profile).
'''
# Modules
import os, sys
import argparse
from collections import defaultdict
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alphal
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.argparse_utils import RoiAction
from hivwholeseq.reference import load_custom_reference, load_custom_alignment



# Functions
def get_alignment_roi(refseq, positions, VERBOSE=0):
    '''Get the alignment corresponding to the coordinates specified
    
    Alignments are much more solid if done gene-by-gene.
    '''
    refname = refseq.name
    start = positions.min()
    end = positions.max() + 1

    # Check all genes
    for fea in refseq.features:
        feaname = fea.id

        if VERBOSE >= 3:
            print 'Trying', feaname,

        # Check for overlap with our roi (assume only one is good)
        # FIXME: deal with exons
        if feaname in ('tat', 'rev'):
            if VERBOSE >= 3:
                print 'not implemented'
            continue

        fea_start = fea.location.nofuzzy_start
        fea_end = fea.location.nofuzzy_end
        # Full coverage
        if (fea_start <= start) and (fea_end >= end):
            ali = load_custom_alignment('HIV1_FLT_2011_'+feaname+'_DNA_subB')

            if VERBOSE >= 3:
                print 'gotcha!'

            break

        if VERBOSE >= 3:
            print 'not good.'


        # TODO: implement partial coverage and missing files

    else:
        raise ValueError('No appropriate genes found for this ROI.')

    # Delete alignment columns that are absent from the chosen reference 
    # the reference must be part of the ali
    for seq in ali:
        if refname in seq.id:
            break
    else:
        raise ValueError('Reference not found in the alignment')

    seqm = np.array(seq, 'S1')
    ind_nongap = (seqm != '-').nonzero()[0]
    mapcoref = np.vstack([fea_start + np.arange(len(ind_nongap)),
                          ind_nongap]).T

    mapcoref = mapcoref[positions - fea_start]

    alim = np.array(ali, 'S1')[:, mapcoref[:, 1]]

    return alim



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Distance btw patient and subtype',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--roi', action=RoiAction, required=True,
                        help='Region of Interest (e.g. F1 100 200)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot results')

    args = parser.parse_args()
    pnames = args.patients
    roi = args.roi
    VERBOSE = args.verbose
    use_plot = args.plot

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]
    
    refname = 'HXB2'
    refseq = load_custom_reference(refname, 'gb')
    refseq.name = refname
    refm = np.array(refseq)

    (fragment, roi_start, roi_end) = roi
    ali = None

    plot_data = defaultdict(list)
    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        if VERBOSE >= 1:
            print pname

        # Find fragment start in reference
        from hivwholeseq.genome_info import find_region_edges
        patseq = patient.get_reference(fragment)
        edges = [str(patseq.seq[:20]), str(patseq.seq[-20:])]
        (frag_start, frag_end) = find_region_edges(refm, edges)


        mapco = patient.get_map_coordinates_reference(fragment, refname=refname,
                                                      roi=roi[1:])

        # Patient afts
        (aft, ind) = patient.get_allele_frequency_trajectories(fragment)
        aft = aft[:, :4, mapco[:, 1]]

        # Subtype B MSA and afs
        ali = get_alignment_roi(refseq,
                                mapco[:, 0],
                                VERBOSE=VERBOSE)
        af = np.vstack([(ali == alpha).mean(axis=0) for alpha in alphal[:4]])

        # FIXME: do better at gaps

        # Calculate distance, e.g. Kullback-Leilbler divergence: it is not symmetric
        # but neither is our data
        divKL = (af * np.log((af + 1e-6) / (aft + 1e-6))).sum(axis=1).sum(axis=1)
        divKLinv = (aft * np.log((aft + 1e-6) / (af + 1e-6))).sum(axis=1).sum(axis=1)

        times = patient.times[ind]

        plot_data['divKL'].append((pname, times, divKL))
        plot_data['divKLinv'].append((pname, times, divKL))

    if use_plot:
        # div KL
        fig, ax = plt.subplots()
        ax.set_xlabel('Time from infection [days]')
        ax.set_ylabel('$D_{KL}($subtype B$ || $patient$)$', fontsize=18)
        
        for i, (pname, x, y) in enumerate(plot_data['divKL']):
            ax.plot(x, y, lw=2, label=pname,
                    color=cm.jet(1.0 * i / len(plot_data['divKL'])))
        
        ax.legend(loc=1, title='Patient:', fontsize=12)
        ax.set_title(' '.join(map(str, roi)))
        ax.grid(True)
        plt.tight_layout()

        # div KL inverse
        fig, ax = plt.subplots()
        ax.set_xlabel('Time from infection [days]')
        ax.set_ylabel('$D_{KL}($patient$ || $subtype B$)$', fontsize=18)
        
        for i, (pname, x, y) in enumerate(plot_data['divKLinv']):
            ax.plot(x, y, lw=2, label=pname,
                    color=cm.jet(1.0 * i / len(plot_data['divKLinv'])))
        
        ax.legend(loc=1, title='Patient:')
        ax.set_title(' '.join(map(str, roi)))
        ax.grid(True)
        plt.tight_layout()

        
        plt.ion()
        plt.show()

