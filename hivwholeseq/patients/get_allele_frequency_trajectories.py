# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/11/13
content:    Collect the allele frequencies of all samples to the initial
            consensus and save them as a matrix into a single trajectory file.
'''
# Modules
import argparse
import numpy as np
from Bio import SeqIO

from hivwholeseq.miseq import alpha
from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file, \
        get_allele_counts_insertions_from_file_unfiltered, \
        filter_nus
from hivwholeseq.patients.patients import get_patient
from hivwholeseq.patients.filenames import get_initial_consensus_filename, \
        get_mapped_to_initial_filename, get_allele_frequency_trajectories_filename
from hivwholeseq.patients.one_site_statistics import plot_allele_frequency_trajectories as plot_nus
from hivwholeseq.patients.one_site_statistics import plot_allele_frequency_trajectories_3d as plot_nus_3d



# Globals



# Functions
def get_allele_frequency_trajectories(patient, fragment, qual_min=35, VERBOSE=0):
    '''Scan the reads of all samples and write to a single file'''
    # Get the initial consensus
    refseq = SeqIO.read(get_initial_consensus_filename(pname, fragment), 'fasta')

    # Prepare output data structures
    nus_traj = np.zeros((len(patient.samples), len(alpha), len(refseq)))
    
    
    for it, sample in enumerate(patient.samples):
        input_filename = get_mapped_to_initial_filename(pname, sample, fragment, type='bam')

        (counts, inserts) = get_allele_counts_insertions_from_file_unfiltered(input_filename,
                                                                   len(refseq),
                                                                   qual_min=qual_min,
                                                                   VERBOSE=VERBOSE)

        # Take the filtered frequencies, blending in the read types
        nu = filter_nus(counts)
        nus_traj[it] = nu

    #FIXME: test, etc.

    return nus_traj    


# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Update initial consensus')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--save', action='store_true',
                        help='Save the allele frequency trajectories to file')
    parser.add_argument('--plot', action='store_true',
                        help='Plot the allele frequency trajectories')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose
    save_to_file = args.save
    plot = args.plot

    # Get patient and the sequenced samples (and sort them by date)
    patient = get_patient(pname)

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Iterate over samples and fragments
    for fragment in fragments:

        # TODO: decide whether or not to move the single allele counts to a
        # separate script, to parallelize (speed up)
            
        nus = np.load(get_allele_frequency_trajectories_filename(pname, fragment))
        #nus = get_allele_frequency_trajectories(patient, fragment, VERBOSE=VERBOSE)

        if save_to_file:
            nus.dump(get_allele_frequency_trajectories_filename(pname, fragment))

        if plot:
            import matplotlib.pyplot as plt

            times = patient.times()
            plot_nus(times, nus, title='Patient '+pname+', '+fragment, VERBOSE=VERBOSE)

            plot_nus_3d(times, nus, title='Patient '+pname+', '+fragment, VERBOSE=VERBOSE)

            plt.tight_layout()
            plt.ion()
            plt.show()
