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
from hivwholeseq.patients.patients import patients as patients_all
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_foldername, get_coordinate_map_filename


# Function
def build_coordinate_map(pname, fragment, refname, VERBOSE=0):
    '''Build the coordinate map'''
    ref_rec = load_custom_reference(refname)
    if 'F5' not in fragment:
        frag_spec = fragment+'o'
    else:
        # Take outermost, F5ao
        frag_spec = 'F5ao'

    if VERBOSE >= 2:
        print frag_spec

    (fragment_start, fragment_end) = find_fragment(ref_rec, frag_spec)
    ref_rec_frag = ref_rec[fragment_start: fragment_end]
    
    ref_init_fn = get_initial_reference_filename(pname, fragment)
    ref_init_rec = SeqIO.read(ref_init_fn, 'fasta')

    from seqanpy import align_global
    ali = align_global(str(ref_rec_frag.seq), str(ref_init_rec.seq))

    # Only identify bijective map
    mapco = []
    pos_ref = fragment_start
    pos_ini = 0
    for col in xrange(len(ali[1])):
        nuc_ref = ali[1][col]
        nuc_ini = ali[2][col]
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
    parser = argparse.ArgumentParser(description='Update initial consensus')
    parser.add_argument('--reference', default='HXB2',
                        help='Select reference strain to align against')
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='+',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--save', action='store_true',
                        help='Save the allele frequency trajectories to file')

    args = parser.parse_args()
    refname = args.reference
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    save_to_file = args.save

    if pnames is None:
        patients = [p for p in patients_all]
    else:
        patients = [p for p in patients_all if p.id in pnames]    
    if VERBOSE >= 3:
        print 'patients', map(attrgetter('id'), patients)
    if not len(patients):
        raise ValueError('No patients found!')

    if fragments is None:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    maps_coord = defaultdict(dict)
    for fragment in fragments:
        for patient in patients:
            pname = patient.id

            mapco = build_coordinate_map(pname, fragment, refname, VERBOSE=VERBOSE)

            maps_coord[fragment][pname] = mapco 

            if save_to_file:
                out_fn = get_coordinate_map_filename(pname, fragment, refname=refname)
                np.savetxt(out_fn, np.array(mapco, int), fmt='%d')
                if VERBOSE:
                    print 'Saved to file:', pname, fragment

