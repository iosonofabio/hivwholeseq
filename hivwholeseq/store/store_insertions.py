#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       26/06/14
content:    Calculate the allele counts and frequencies for each patient sample
            (and both PCR1 and PCR2 if present).
'''
# Modules
import os
from warnings import warn
import argparse
import numpy as np
from Bio import SeqIO

from hivwholeseq.utils.argparse import PatientsAction
from hivwholeseq.utils.exceptions import NoDataWarning
from hivwholeseq.patients.samples import load_samples_sequenced as lssp
from hivwholeseq.patients.samples import SamplePat
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_mapped_filtered_filename, get_insertions_filename
from hivwholeseq.utils.one_site_statistics import get_allele_counts_insertions_from_file as gac
from hivwholeseq.cluster.fork_cluster import fork_get_insertions_patient as fork_self 



# Functions
def save_insertions(filename, insertions):
    '''Save insertions to file'''
    import cPickle as pickle
    with open(filename, 'w') as f:
        pickle.dump(insertions, f, protocol=-1)



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Store insertions',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=True)
    pats_or_samples.add_argument('--patients', action=PatientsAction,
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to analyze')
    parser.add_argument('--fragments', nargs='+',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=2,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save the allele counts to file')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--qualmin', type=int, default=30,
                        help='Minimal quality of base to call')
    parser.add_argument('--PCR', type=int, default=1,
                        help='Analyze only reads from this PCR (e.g. 1)')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    fragments = args.fragments
    submit = args.submit
    VERBOSE = args.verbose
    save_to_file = args.save
    qual_min = args.qualmin
    PCR = args.PCR

    samples = lssp()
    if pnames is not None:
        samples = samples.loc[samples.patient.isin(pnames)]
    elif samplenames is not None:
        samples = samples.loc[samples.index.isin(samplenames)]

    if VERBOSE >= 2:
        print 'samples', samples.index.tolist()

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    for fragment in fragments:
        inses = []
        for samplename, sample in samples.iterrows():
            if submit:
                fork_self(samplename, fragment, VERBOSE=VERBOSE, qual_min=qual_min)
                continue

            if VERBOSE >= 1:
                print fragment, samplename

            sample = SamplePat(sample)
            pname = sample.patient
            refseq = SeqIO.read(get_initial_reference_filename(pname, fragment), 'fasta')

            fn = sample.get_mapped_filtered_filename(fragment, PCR=PCR)
            if not os.path.isfile(fn):
                warn('No BAM file found', NoDataWarning)
                continue

            _, inse = gac(fn, len(refseq), qual_min=qual_min, VERBOSE=VERBOSE)
            inses.append(inse)

            if save_to_file:
                fn_out = sample.get_insertions_filename(fragment, PCR=PCR,
                                                        qual_min=qual_min)
                save_insertions(fn_out, inse)

                if VERBOSE >= 2:
                    print 'Insertions saved:', samplename, fragment
