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

from hivwholeseq.sequencing.primer_info import find_fragment
from hivwholeseq.reference import load_custom_reference
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_foldername, get_coordinate_map_filename


# Function
def build_coordinate_map(refseq, patseq, VERBOSE=0):
    '''Build the coordinate map'''
    from seqanpy import align_overlap
    (score, ali1, ali2) = align_overlap(refseq, patseq)
    patseq_start = len(ali2) - len(ali2.lstrip('-'))
    patseq_end = len(ali2.rstrip('-'))

    if VERBOSE >= 3:
        from hivwholeseq.sequence_utils import pretty_print_pairwise_ali
        pretty_print_pairwise_ali([ali1[patseq_start: patseq_end],
                                   ali2[patseq_start: patseq_end]],
                                  name1=refseq.name, name2=patseq.name)

    # Only identify bijective map
    mapco = []
    pos_ref = patseq_start
    pos_ini = 0
    for col in xrange(patseq_start, patseq_end):
        nuc_ref = ali1[col]
        nuc_ini = ali2[col]
        if (nuc_ref != '-') and (nuc_ini != '-'):
            mapco.append((pos_ref, pos_ini))
            pos_ref += 1
            pos_ini += 1
        elif (nuc_ref != '-'):
            pos_ref += 1
        elif (nuc_ini != '-'):
            pos_ini += 1

    return mapco



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Map patient coordinates to reference',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--reference', default='HXB2',
                        help='Select reference strain to align against')
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='+',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--save', action='store_true',
                        help='Save the map to file')

    args = parser.parse_args()
    refname = args.reference
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    save_to_file = args.save

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]
    if VERBOSE >= 3:
        print 'patients', patients.index
    if not len(patients):
        raise ValueError('No patients found!')

    if fragments is None:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    maps_coord = defaultdict(dict)
    for fragment in fragments:
        refseq = load_custom_reference(refname)
        for pname, patient in patients.iterrows():
            if VERBOSE >= 1:
                print fragment, pname
            patient = Patient(patient)
            patseq = patient.get_reference(fragment)

            mapco = build_coordinate_map(refseq, patseq, VERBOSE=VERBOSE)

            maps_coord[fragment][pname] = mapco 

            if save_to_file:
                out_fn = get_coordinate_map_filename(pname, fragment, refname=refname)
                np.savetxt(out_fn, np.array(mapco, int), fmt='%d')
                if VERBOSE:
                    print 'Saved to file:', pname, fragment

