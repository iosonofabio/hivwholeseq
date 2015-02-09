# vim: fdm=marker
'''
author:     Fabio Zanini
date:       05/02/15
content:    Try to date haplotypes from cell samples in terms of previous
            RNA samples.
'''
# Modules
import os, sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

from hivwholeseq.utils.argparse import RoiAction
from hivwholeseq.patients.patients import load_patient
from hivwholeseq.patients.samples import SamplePat



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Date the cell haplotypes',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patient', required=True,
                        help='Patient to study')
    parser.add_argument('--roi', required=True, action=RoiAction,
                        help='Region of interest (e.g. F1 300 350)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    pname = args.patient
    roi = args.roi
    VERBOSE = args.verbose

    if VERBOSE >= 1:
        print pname, roi

    patient = load_patient(pname)
    roi = patient.get_fragmented_roi(roi, include_genomewide=False)

    hct, ind, seqs = patient.get_local_haplotype_count_trajectories(*roi, VERBOSE=VERBOSE)

    if patient.samples.iloc[ind[-1]]['sample type'] != 'PBMC':
        raise ValueError('Cell sample not found')

    hindcell = hct[-1].argsort()[::-1][:(hct[-1] > 1).sum()]
    hft = (1.0 * hct.T / hct.sum(axis=1)).T
    hind0 = hct[0].argmax()

    print 'Patient:', pname
    print 'Region:', roi
    print 'Initial consensus:'
    print seqs[hind0]
    print 
    print 'Rank in cell                   time [days from infection]'
    print '\t'.join(['            '] + ['{:>4d}'.format(int(t)) for t in patient.times[ind[:-1]]] + ['cell'])
    print '-' * 110
    for ii, i in enumerate(hindcell, 1):
        print 'Rank', '{:>2d}'.format(ii), '\t', '\t'.join(map('{:>4.1%}'.format, hft[:, i]))
