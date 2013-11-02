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

from mapping.miseq import alpha
from mapping.one_site_statistics import get_allele_counts_insertions_from_file, \
        filter_nus
from mapping.patients.patients import get_patient, get_sequenced_samples
from mapping.patients.filenames import get_initial_consensus_filename, \
        get_mapped_to_initial_filename



# Globals



# Functions
def get_allele_frequency_trajectories(patient, fragment, qual_min=35, VERBOSE=0):
    '''Scan the reads of all samples and write to a single file'''
    # Get the initial consensus
    refseq = SeqIO.read(get_initial_consensus_filename(pname, fragment), 'fasta')

    # Prepare output data structures
    nus_traj = np.zeros((len(patient['samples']), len(alpha), len(refseq)))
    
    
    for it, sample in enumerate(patient['samples']):
        # The mapped reads should be already clean at this point
        input_filename = get_mapped_to_initial_filename(pname, sample, fragment,
                                                         type='bam')
        (counts, inserts) = get_allele_counts_insertions_from_file(input_filename,
                                                                   len(refseq),
                                                                   qual_min=qual_min,
                                                                   VERBOSE=VERBOSE)

        # Take the filtered frequencies, blending in the read types
        nu = filter_nus(counts, counts.sum(axis=0))
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
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose

    # Get patient and the sequenced samples (and sort them by date)
    patient = get_patient(pname)
    patient['samples'] = get_sequenced_samples(patient)

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Iterate over samples and fragments
    for fragment in fragments:

        # Fragment F5 has two sets of primers, so it's a mess in this case,
        # because the initial consensus might not have all the bases :-/
        # FIXME: make sure the initial consensus has ;-)
        if fragment == 'F5':
            raise ValueError('F5a/F5b harmony not implemented yet!')

        # TODO: decide whether or not to move the single allele counts to a
        # separate script, to parallelize (speed up)
            
        get_allele_frequency_trajectories(patient, fragment, VERBOSE=VERBOSE)
