# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/01/15
content:    Build an alignment of all patient references, to check they are all
            fine, starting and ending from the same coordinates, and so on.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import AlignIO

from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.utils.sequence import align_muscle
from hivwholeseq.patients.filenames import get_reference_alignment_filename


# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Align references',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--fragments', nargs='+',
                        default=['F'+str(i) for i in xrange(1, 7)] + ['genomewide'],
                        help='Fragments to taylor')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    VERBOSE = args.verbose
    fragments = args.fragments

    patients = load_patients()

    for fragment in fragments:
        if VERBOSE >= 1:
            print fragment

        refs = []
        for pname, patient in patients.iterrows():
            if VERBOSE >= 2:
                print pname
            patient = Patient(patient)
            refs.append(patient.get_reference(fragment))

        ali = align_muscle(*refs, sort=True)

        # Check whether all references are complete (using the longest primers)
        if VERBOSE >= 2:
            print 'Check alignment'
        alim = np.array(ali)
        if (alim[:, :4] == '-').any():
            raise ValueError('Gaps at the beginning of fragment found')
        elif (alim[:, -4:] == '-').any():
            raise ValueError('Gaps at the end of fragment found')

        if VERBOSE >= 2:
            print 'Save to file'
        fn = get_reference_alignment_filename(fragment)
        AlignIO.write(ali, fn, 'fasta')

 
