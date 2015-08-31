# vim: fdm=marker
'''
author:     Fabio Zanini
date:       31/08/15
content:    Store local divergence and diversity in a sliding window.
'''
# Modules
import os
import sys
import argparse
from operator import itemgetter
import numpy as np
from Bio import SeqIO

from hivwholeseq.miseq import alpha
from hivwholeseq.utils.argparse import PatientsAction
from hivwholeseq.patients.patients import load_patients, Patient


# Functions
def get_divergence_diversity_sliding(aft, block_length, VERBOSE=0):
    '''Get local divergence and diversity in a sliding window'''
    cons_ind = Patient.get_initial_consensus_noinsertions(aft, return_ind=True)
    ind_N = cons_ind == 5
    cons_ind[ind_N] = 0
    aft_nonanc = 1.0 - aft[:, cons_ind, np.arange(aft.shape[2])]
    aft_nonanc[:, ind_N] = 0

    aft_var = (aft * (1 - aft)).sum(axis=1)

    struct = np.ones(block_length)

    dg = np.ma.array(np.apply_along_axis(lambda x: np.convolve(x, struct, mode='valid'),
                                         axis=1, arr=aft_nonanc), hard_mask=True)
    ds = np.ma.array(np.apply_along_axis(lambda x: np.convolve(x, struct, mode='valid'),
                                         axis=1, arr=aft_var), hard_mask=True)

    # NOTE: normalization happens based on actual coverage
    norm = np.apply_along_axis(lambda x: np.convolve(x, struct, mode='valid'),
                               axis=1, arr=(-aft[:, 0].mask))

    dg.mask = norm < block_length
    dg /= norm

    ds.mask = norm < block_length
    ds /= norm

    x = np.arange(dg.shape[1]) + (block_length - 1) / 2.0

    return (x, dg, ds)


def get_divergence_diversity_blocks(aft, block_length, VERBOSE=0):
    '''Get local divergence and diversity in blocks'''
    cons_ind = aft[0].argmax(axis=0)
    n_blocks = aft.shape[2] // block_length
    dg = np.zeros((len(aft), n_blocks))
    ds = np.zeros_like(dg)
    for n_block in xrange(n_blocks):
        for pos in xrange(block_length):
            pos += n_block * block_length
            af = aft[:, :, pos]
            af_nonanc = 1.0 - af[:, cons_ind[pos]]
            dg[:, n_block] += af_nonanc
            ds[:, n_block] += (af * (1 - af)).sum(axis=1)
    dg /= block_length
    ds /= block_length

    dg = np.ma.array(dg)
    ds = np.ma.array(ds)

    x = (np.arange(n_blocks) + 0.5) * block_length

    return (x, dg, ds)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Store local divergence and diversity',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', action=PatientsAction,
                        help='Patients to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=2,
                        help='Verbosity level [0-4]')
    parser.add_argument('--block-length', type=int, default=150,
                        help='Length of block to consider')
    parser.add_argument('--sliding', action='store_true',
                        help='Use a sliding window instead of dividing into blocks')
    parser.add_argument('--save', action='store_true',
                        help='Save the divergence and diversity to file')

    args = parser.parse_args()
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    block_length = args.block_length
    use_sliding = args.sliding
    save_to_file = args.save

    patients = load_patients()
    if pnames is not None:
        patients = patients.iloc[patients.index.isin(pnames)]

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    dgs = {}
    dss = {}

    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()

        for ifr, fragment in enumerate(fragments):
            if VERBOSE >= 1:
                print pname, fragment

            aft, ind = patient.get_allele_frequency_trajectories(fragment,
                                                                 cov_min=100)

            # NOTE: Ns should be excluded from diversity and divergence
            aft = aft[:, :5, :]

            if use_sliding:
                (x, dg, ds) = get_divergence_diversity_sliding(aft, block_length,
                                                               VERBOSE=VERBOSE)
            else:
                (x, dg, ds) = get_divergence_diversity_blocks(aft, block_length,
                                                              VERBOSE=VERBOSE)

            # FIXME: avoid this var to get different conv and aft indices
            times = patient.times[ind]
            dgs[(pname, fragment)] = (patient.times[ind], dg)
            dss[(pname, fragment)] = (patient.times[ind], ds)

            if save_to_file:
                from hivwholeseq.patients.filenames import \
                        get_divergence_trajectories_local_filename, \
                        get_diversity_trajectories_local_filename

                # savez does not support masked arrays
                dg_save = np.array(dg).copy()
                dg_save[dg.mask] = -1
                ds_save = np.array(ds).copy()
                ds_save[ds.mask] = -1
                L = aft.shape[2]

                fn_out = get_divergence_trajectories_local_filename(pname, fragment)
                np.savez(fn_out, ind=ind, dg=dg_save, L=L, block_length=[block_length])
                fn_out = get_diversity_trajectories_local_filename(pname, fragment)
                np.savez(fn_out, ind=ind, ds=ds_save, L=L, block_length=[block_length])
                if VERBOSE >= 1:
                    print 'saved to file'
