# vim: fdm=marker
'''
author:     Fabio Zanini
date:       22/05/14
content:    Build a coordinate map of the initial reference of a patient to an
            external reference seq (e.g. HXB2). This is useful to quickly find
            genes and stuff like that.
'''
# Modules
import os
import argparse
from operator import attrgetter
from collections import defaultdict
import numpy as np
from Bio import SeqIO

from hivwholeseq.reference import load_custom_reference
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_foldername, get_coordinate_map_filename



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Map patient coordinates to reference',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--reference', default='HXB2',
                        help='Select reference strain to align against')
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--regions', nargs='+',
                        help='regions to make coordinate maps for (e.g. V3 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--save', action='store_true',
                        help='Save the map to file')

    args = parser.parse_args()
    refname = args.reference
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    save_to_file = args.save

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]
    if VERBOSE >= 3:
        print 'patients', patients.index
    if not len(patients):
        raise ValueError('No patients found!')

    maps_coord = defaultdict(dict)
    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        # Make maps for all annotations if not explicit
        if regions is None:
            patseqann = patient.get_reference('genomewide', format='gb')
            regionspat = map(attrgetter('id'), patseqann.features) + ['genomewide']
        else:
            regionspat = regions

        for region in regionspat:
            if VERBOSE >= 1:
                print pname, region

            coomap = patient.get_map_coordinates_reference(region,
                                                           refname=refname)

