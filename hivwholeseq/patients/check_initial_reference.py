# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/06/14
content:    Check whether we need a new initial reference. Reasons might be a new
            sample sequenced that comes before all current ones, or the current
            reference has some genes not properly assembled.
'''
# Modules
import os
import argparse

from hivwholeseq.patients.patients import load_patients, load_patient, Patient
from hivwholeseq.patients.filenames import get_sample_foldername


# Functions
def check_sample_folders(pname, samplename_init, VERBOSE=0):
    '''Check presence of initial sample among the sample folders'''
    for PCR in (1, 2):
        subdirname = get_sample_foldername(pname, samplename_init+'_PCR'+str(PCR))
        if VERBOSE >= 2:
            print pname, 'trying', subdirname
        if os.path.isdir(subdirname):
            return True
    return False




# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Check patient samples')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')

    args = parser.parse_args()
    VERBOSE = args.verbose
    pnames = args.patients

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[patients.index.isin(pnames)]

    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        if VERBOSE:
            print patient.name

        patient.discard_nonsequenced_samples()

        # Check whether patient has at least three time points, else ignore
        if patient.samples.shape[0] < 3:
            print 'WARNING: patient has less than three samples sequenced. Skipping.'
            continue
    
        # Check whether the first sample has a folder already. If not, it's
        # not going to work
        samplename_init = str(patient.initial_sample.name)
        check = check_sample_folders(patient.name, samplename_init, VERBOSE=VERBOSE)
        if not check:
            print 'ERROR: Folder for initial sample, '+samplename_init+', not found.'
        elif VERBOSE:
            print 'OK: Folder for initial sample, '+samplename_init+', found.'

