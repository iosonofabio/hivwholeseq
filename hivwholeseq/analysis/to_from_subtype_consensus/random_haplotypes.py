# vim: fdm=marker
'''
author:     Fabio Zanini
date:       12/02/15
content:    Take random haplotypes at different time points and see how many
            sites they have that differ from the subtype consensus.
'''
# Modules
import os, sys
import argparse
from collections import defaultdict
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alphal, alpha
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.reference import load_custom_reference, load_custom_alignment
from hivwholeseq.cross_sectional.get_subtype_consensus import (
    get_subtype_reference_alignment_consensus)


# Globals
refname = 'HXB2'
Sbins = [0, 0.05, 0.2, 0.5, 2]



# Functions
def align_to_reference(conss, alism, VERBOSE=0):
    '''Align alignment to reference'''
    from hivwholeseq.utils.sequence import align_muscle

    # The consensus cannot have gaps
    conss = conss.replace('-', 'N')

    alim = np.array(align_muscle(*(np.append(conss, alism)), sort=True), 'S1')
    alism = np.array(map(''.join, alim[1:, alim[0] != '-']))
    return alism



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
    

    data = []
    for region in regions:
        if VERBOSE >= 1:
            print region

        # Load subtype alignment consensus
        consrec = get_subtype_reference_alignment_consensus(region)
        consm = np.array(consrec)
        conss = ''.join(consrec)

        for pname, patient in patients.iterrows():
            if VERBOSE >= 1:
                print pname, region
            patient = Patient(patient)

            if VERBOSE >= 2:
                print 'Get haplotype count'
            hct, ind, alim = patient.get_haplotype_count_trajectory(region)
            alism = np.array(map(''.join, alim), 'S'+str(len(alim[0])))

            # Exclude time points without counts
            ii = (hct > 0).any(axis=1)
            hct = hct[ii]
            ind = ind[ii]

            if VERBOSE >= 2:
                print 'Align to reference'
            alisma = align_to_reference(conss, alism, VERBOSE=VERBOSE)
            alima = np.array([np.fromstring(seq, 'S1') for seq in alisma])

            if VERBOSE >= 2:
                print 'Get frequencies'
            hft = (1.0 * hct.T / hct.sum(axis=1)).T

            if VERBOSE >= 2:
                print 'Get initial consensus'
            consm0 = alima[hft[0].argmax()]

            if VERBOSE >= 2:
                print 'Get cumulative frequencies'
            hcft = hft.cumsum(axis=1)

            if VERBOSE >= 2:
                print 'Get random haplotypes'
            rs = np.random.rand(ind.shape[0], 1000)
            haplors = []
            for it, rt in enumerate(rs):
                for j, r in enumerate(rt):
                    jseq = (r <= hcft[it]).nonzero()[0][0]
                    haplors.append(alima[jseq])
            haplors = np.array(haplors).reshape(list(rs.shape)+[alima.shape[1]])

            if VERBOSE >= 2:
                print 'Get distributions'
            frac_consa = (haplors == consm).mean(axis=-1)
            frac_consa_sort = np.array(map(np.sort, frac_consa))

            frac_cons0 = (haplors == consm0).mean(axis=-1)
            frac_cons0_sort = np.array(map(np.sort, frac_cons0))

            if use_plot:
                if VERBOSE >= 2:
                    print 'Plot'

                fig, ax = plt.subplots()

                for ii, i in enumerate(ind):
                    time = patient.times[i]

                    # Plot distribution from subtype consensus
                    da = frac_consa_sort[ii]
                    ax.plot(da, 1 - np.linspace(0, 1, len(da)),
                            ls='-',
                            lw=2,
                            color=cm.jet(1.0 * i / ind.max()))

                    # Plot distribution from initial consensus
                    d0 = frac_cons0_sort[ii]
                    ax.plot(d0, 1 - np.linspace(0, 1, len(d0)),
                            ls='--',
                            lw=2,
                            color=cm.jet(1.0 * i / ind.max()))

                ax.set_xlabel('Fraction of matching sites')
                ax.set_ylabel('Fraction of variants > x')
                ax.set_xlim(0.85, 1.01)
                ax.set_ylim(-0.05, 1.05)
                ax.grid(True)

                plt.tight_layout()
                plt.ion()
                plt.show()
