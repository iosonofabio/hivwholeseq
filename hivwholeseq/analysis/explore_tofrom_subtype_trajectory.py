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
from hivwholeseq.utils.argparse import RoiAction
from hivwholeseq.reference import load_custom_reference, load_custom_alignment
from hivwholeseq.cross_sectional.get_subtype_entropy import (
    get_subtype_reference_alignment)


# Globals
refname = 'HXB2'



# Functions



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Mutations away from/towards subtype',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Genomic regions (e.g. V3 IN)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot results')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    use_plot = args.plot

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]
    

    for region in regions:

        # Load HXB2
        refseq = load_custom_reference(refname, 'gb', region=region)
        refseq.name = refname
        refm = np.array(refseq)

        # Load subtype alignment
        ali = get_subtype_reference_alignment(region)

    # FIXME
    sys.exit()


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

