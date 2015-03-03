# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/09/14
content:    Grab all consensi from a patient and align them.
'''
# Modules
import sys
import os
import argparse
import numpy as np

from hivwholeseq.patients.patients import load_patients, Patient



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Align consensi',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--regions', nargs='*',
                        help='Regions to analyze (e.g. V3 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose

    patients = load_patients()
    if pnames:
        patients = patients.loc[pnames]

    alis = {}
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()

        # Guess regions if not specified
        if regions is None:
            refseq_gw = patient.get_reference('genomewide', 'gb')
            regionspat = map(attrgetter('id'), refseq_gw.features) + ['genomewide']
        else:
            regionspat = regions

        for region in regionspat:
            if VERBOSE >= 1:
                print pname, region

            ali = patient.get_consensi_alignment(region)

            alis[(region, pname)] = ali
