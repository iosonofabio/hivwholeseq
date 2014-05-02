# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/03/14
content:    Merge allele frequency trajectories of all fragments.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO

from seqanpy import align_overlap

from hivwholeseq.miseq import alpha
from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file, \
        get_allele_counts_insertions_from_file_unfiltered, \
        filter_nus
from hivwholeseq.patients.patients import get_patient
from hivwholeseq.patients.filenames import get_initial_consensus_filename, \
        get_mapped_to_initial_filename, get_allele_frequency_trajectories_filename, \
        get_allele_count_trajectories_filename
from hivwholeseq.patients.one_site_statistics import plot_allele_frequency_trajectories as plot_nus
from hivwholeseq.patients.one_site_statistics import plot_allele_frequency_trajectories_3d as plot_nus_3d



# Functions
def merge_allele_frequency_trajectories(conss_genomewide, acss, VERBOSE=0):
    '''Merge the allele frequency trajectories of all fragments by summing counts'''

    consm = np.fromstring(conss_genomewide, 'S1')
    acs = np.zeros((acss[0].shape[0], len(alpha), len(conss_genomewide)), int)

    # Put the first fragment
    pos_ref = 0
    for pos in xrange(acss[0].shape[2]):
        acs[:, :, pos] = acss[0][:, :, pos]

    # Blend in subsequent fragments
    for ifr, acsi in enumerate(acss[1:], 1):
        seed = alpha[acsi[0, :, :30].argmax(axis=0)]
        sl = len(seed)
        n_match = np.array([(consm[i: i + sl] == seed).sum()
                            for i in xrange(pos_ref, len(consm) - sl)], int)
        pos_seed = n_match.argmax()
        if n_match[pos_seed] < 0.75 * sl:
            raise ValueError('Fragment could not be found')
        pos_seed += pos_ref

        if VERBOSE >= 2:
            print 'F'+str(ifr+1), 'seed:', pos_seed

        # Align the afc block to the global consensus
        cons_acsi = ''.join(alpha[acsi[0].argmax(axis=0)])
        ali = align_overlap(conss_genomewide[pos_seed: pos_seed + acsi.shape[2] + 100],
                            cons_acsi, band=100)
        ali = list(ali)
        ali[2] = ali[2].rstrip('-')
        ali[1] = ali[1][:len(ali[2])]
        # Scan the new fragment and check the major allele, infer gaps that way
        # (somewhat sloppy, but ok)
        pos_gw = pos_seed
        pos_acsi = 0
        for pos in xrange(len(ali[2])):
            # Match
            if (ali[1][pos] != '-') and (ali[2][pos] != '-'):
                acs[:, :, pos_gw] += acsi[:, :, pos_acsi]
                pos_gw += 1
                pos_acsi += 1
            
            # Insert in the block: discard
            elif (ali[1][pos] == '-'):
                pos_acsi += 1

            # Deletion in block: ignore (we should dig into inserts, but as far
            # as frequencies are concerned, it should be fine)
            else:
                pos_gw += 1

        pos_ref = pos_seed

    # Normalize to frequencies
    afs = (1.0 * acs.swapaxes(0, 1) / acs.sum(axis=1)).swapaxes(0, 1)

    return afs



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Update initial consensus')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--save', action='store_true',
                        help='Save the allele frequency trajectories to file')
    parser.add_argument('--plot', action='store_true',
                        help='Plot the allele frequency trajectories')

    args = parser.parse_args()
    pname = args.patient
    VERBOSE = args.verbose
    save_to_file = args.save
    plot = args.plot

    # Get patient and the sequenced samples (and sort them by date)
    patient = get_patient(pname)

    # If the script is called with no fragment, iterate over all
    fragments = ['F'+str(i) for i in xrange(1, 7)]
    
    # Collect the allele count trajectories
    acss = []
    for fragment in fragments:
        act_filename = get_allele_count_trajectories_filename(pname, fragment)
        acsi = np.load(act_filename)
        acss.append(acsi)

    # Get the initial genomewide consensus
    cons_filename = get_initial_consensus_filename(pname, 'genomewide')
    cons_rec = SeqIO.read(cons_filename, 'fasta')
    conss_genomewide = str(cons_rec.seq)

    # Merge allele counts
    afs = merge_allele_frequency_trajectories(conss_genomewide, acss, VERBOSE=VERBOSE)

    if save_to_file:
        afs.dump(get_allele_frequency_trajectories_filename(pname, 'genomewide'))

    if plot:
        import matplotlib.pyplot as plt
        times = patient.times()
        plot_nus(times, afs, title='Patient '+pname, VERBOSE=VERBOSE)
        plot_nus_3d(times, afs, title='Patient '+pname, VERBOSE=VERBOSE)

        plt.ion()
        plt.show()
    
    ## Average in windows of 50 bp
    #window_size = 50
    #afs_win = np.zeros((afs.shape[0], afs.shape[1], afs.shape[2] // window_size))
    #for i in xrange(afs_win.shape[2]):
    #    afs_win[:, :, i] = afs[:, :, i * window_size: (i+1) * window_size].mean(axis=2)

    #if plot:
    #    times = patient.times()
    #    plot_nus_3d(times, afs_win, title='Patient '+pname, VERBOSE=VERBOSE)

