#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       03/12/14
content:    Compute alignments of haplotypes in a few regions and store them for
            the website.
'''
# Modules
import os
import sys
import argparse
import numpy as np
import StringIO
from Bio import AlignIO

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.sequence_utils import build_msa_haplotypes as build_msa
from hivwholeseq.website.filenames import get_precompiled_alignments_filename
from hivwholeseq.fork_cluster import fork_store_haplotypes_website as fork_self


# Globals
regions_all = ['PR', 'V3', 'psi']



# Functions
def store_alignments(alis, filename):
    '''Store alignments to file'''
    import zipfile, zlib
    with zipfile.ZipFile(filename, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
        for i, ali in enumerate(alis):
            f = StringIO.StringIO()
            AlignIO.write(ali['ali'], f, 'fasta')
            zf.writestr(str(ali['time'])+'_days.fasta', f.getvalue())



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Perform PCA on the data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--regions', nargs='+', default=regions_all,
                        help='Regions to store (e.g. V3 PR)')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    submit = args.submit

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    if submit:
        for pname, patient in patients.iterrows():
            for region in regions:
                fork_self(pname, region, VERBOSE=1)
        sys.exit()

    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        for region in regions:
            print patient.code, patient.name, region

            # Compute alignments
            print 'Get haplotypes'
            haplos = patient.get_local_haplotype_trajectories(region, 0, '+oo')
            print 'Align'
            alis = [{'time': patient.times[it],
                     'ali': build_msa(h)}
                    for h, it in zip(*haplos)]

            # Write to file
            print 'Write output file'
            fn_out = get_precompiled_alignments_filename(patient.code, region)
            mkdirs(os.path.dirname(fn_out))
            store_alignments(alis, fn_out)


