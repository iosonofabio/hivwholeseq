#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/12/14
content:    Store local haplotypes.
'''
# Modules
import os
import argparse
from operator import itemgetter, attrgetter
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from Bio import Phylo

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.utils.mapping import align_muscle
from hivwholeseq.patients.patients import load_patients, Patient



# Functions



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Store local haplotypes',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Genomic regions (e.g. V3 IN)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of reads analyzed per sample')
    parser.add_argument('--freqmin', type=int, default=0.01,
                        help='Minimal frequency to keep the haplotype')
    parser.add_argument('--countmin', type=int, default=3,
                        help='Minimal observations to keep the haplotype')
    parser.add_argument('--save', action='store_true',
                        help='Save alignment to file')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    maxreads = args.maxreads
    freqmin = args.freqmin
    countmin = args.countmin
    use_save = args.save

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        for region in regions:
            if VERBOSE >= 1:
                print pname, region
    
            if VERBOSE >= 2:
                print 'Get region haplotypes'
            filters = ['noN', 'mincount='+str(countmin)]
            (hct, ind, seqs) = \
                    patient.get_local_haplotype_count_trajectories(region,
                                                                   filters=filters,
                                                                   VERBOSE=VERBOSE)

            # Filter out tiny-frequency variants
            hft = (1.0 * hct.T / hct.sum(axis=1)).T
            ind_seqs = (hft > freqmin).any(axis=0)
            hft = hft[:, ind_seqs]
            seqs = seqs[ind_seqs]

            # Align sequences
            ali = np.array(align_muscle(*seqs, sort=True))

            if use_save:
                if VERBOSE >= 2:
                    print 'Save to file'
                fn_out = patient.get_haplotype_count_trajectory_filename(region)
                mkdirs(os.path.dirname(fn_out))
                np.savez_compressed(fn_out,
                                    hct=hct,
                                    ind=ind,
                                    times=patient.times[ind],
                                    seqs=seqs,
                                    ali=ali)
