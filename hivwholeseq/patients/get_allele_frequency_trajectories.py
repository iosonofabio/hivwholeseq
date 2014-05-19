#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/11/13
content:    Collect the allele frequencies of all samples to the initial
            consensus and save them as a matrix into a single trajectory file.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO

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
from hivwholeseq.fork_cluster import fork_get_allele_frequency_trajectory as fork_self



# Functions
def get_allele_frequency_trajectories(patient, fragment, qual_min=35, VERBOSE=0):
    '''Scan the reads of all samples and write to a single file'''
    # Get the initial consensus
    refseq = SeqIO.read(get_initial_consensus_filename(pname, fragment), 'fasta')

    # Prepare output data structures
    cos_traj = np.zeros((len(patient.samples), len(alpha), len(refseq)), int)
    nus_traj = np.zeros((len(patient.samples), len(alpha), len(refseq)))
    
    
    for it, sample in enumerate(patient.samples):
        input_filename = get_mapped_to_initial_filename(pname, sample, fragment, type='bam')

        (counts, inserts) = get_allele_counts_insertions_from_file_unfiltered(input_filename,
                                                                   len(refseq),
                                                                   qual_min=qual_min,
                                                                   VERBOSE=VERBOSE)
        # Take the total counts, blending in the read types
        cou = counts.sum(axis=0)
        cos_traj[it] = cou

        # Take the filtered frequencies, blending in the read types
        nu = filter_nus(counts)
        nus_traj[it] = nu

    #FIXME: test, etc.

    return (cos_traj, nus_traj)


# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get allele frequency trajectories')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--save', action='store_true',
                        help='Save the allele frequency trajectories to file')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the allele frequency trajectories')
    parser.add_argument('--logit',nargs='?', const=True, default=False,
                        help='logit scale (log(x/(1-x))')
    parser.add_argument('--PCR1', action='store_true',
                        help='Show only PCR1 samples where possible (still computes all)')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose
    save_to_file = args.save
    plot = args.plot
    submit = args.submit
    use_PCR1 = args.PCR1

    patient = get_patient(pname)
    times = patient.times()
    samplenames = patient.samples
    if use_PCR1:
        # Keep PCR2 only if PCR1 is absent
        ind = np.nonzero(map(lambda x: ('PCR1' in x[1]) or ((times == times[x[0]]).sum() == 1), enumerate(samplenames)))[0]
        times = times[ind]

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Iterate over samples and fragments
    for fragment in fragments:

        # Submit to the cluster self if requested (--save is assumed)
        if submit:
            fork_self(pname, fragment,
                      VERBOSE=VERBOSE)
            continue

        if VERBOSE >= 1:
            print fragment

        # Save new computation?
        act_filename = get_allele_count_trajectories_filename(pname, fragment)
        aft_filename = get_allele_frequency_trajectories_filename(pname, fragment)
        if save_to_file or (not os.path.isfile(aft_filename)):
            (act, aft) = get_allele_frequency_trajectories(patient, fragment, VERBOSE=VERBOSE)
        else:
            aft = np.load(aft_filename)
            act = np.load(act_filename)

        aft[np.isnan(aft)] = 0
        aft[(aft < 1e-5) | (aft > 1)] = 0

        if save_to_file:
            act.dump(get_allele_count_trajectories_filename(pname, fragment))
            aft.dump(get_allele_frequency_trajectories_filename(pname, fragment))

        if use_PCR1:
            aft = aft[ind]
            act = act[ind]

        if plot is not None:
            import matplotlib.pyplot as plt

            if plot in ('2D', '2d', ''):
                plot_nus(times, aft, title='Patient '+pname+', '+fragment, VERBOSE=VERBOSE, logit = args.logit)
            elif plot in ('3D', '3d'):
                plot_nus_3d(times, aft, title='Patient '+pname+', '+fragment, VERBOSE=VERBOSE, logit = args.logit)

    if plot is not None:
        plt.tight_layout()
        plt.ion()
        plt.show()
