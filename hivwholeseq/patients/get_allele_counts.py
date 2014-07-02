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
import argparse
import numpy as np
from Bio import SeqIO

from hivwholeseq.patients.filenames import get_initial_consensus_filename, \
        get_mapped_filtered_filename, get_allele_counts_filename
from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file as gac
from hivwholeseq.patients.patients import load_samples_sequenced as lssp
from hivwholeseq.fork_cluster import fork_get_allele_counts_patient as fork_self 



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get allele frequency trajectories',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=True)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to map (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save the allele frequency trajectories to file')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the allele frequency trajectories')
    parser.add_argument('--logit', action='store_true',
                        help='use logit scale (log(x/(1-x)) in the plots')
    parser.add_argument('--qualmin', type=int, default=30,
                        help='Minimal quality of base to call')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    fragments = args.fragments
    submit = args.submit
    VERBOSE = args.verbose
    save_to_file = args.save
    qual_min = args.qualmin

    samples_pat = lssp()
    if pnames is not None:
        samples_pat = samples_pat.loc[samples_pat.patient.isin(pnames)]
    elif samplenames is not None:
        samples_pat = samples_pat.loc[samples_pat.index.isin(samplenames)]

    if VERBOSE >= 2:
        print 'samples', samples_pat.index.tolist()

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    counts_all = []
    for fragment in fragments:
        counts = []
        for samplename_pat, sample_pat in samples_pat.iterrows():
            if submit:
                fork_self(samplename_pat, fragment, VERBOSE=VERBOSE, qual_min=qual_min)
                continue

            pname = sample_pat.patient
            refseq = SeqIO.read(get_initial_consensus_filename(pname, fragment), 'fasta')
            counts_sample = [None, None]

            # Look for both PCR1 and PCR2, and see what you find
            for PCR in (1, 2):
                fn_out = get_allele_counts_filename(pname, samplename_pat,
                                                    fragment, PCR=PCR, qual_min=qual_min)
                fn = get_mapped_filtered_filename(pname, samplename_pat,
                                                  fragment, PCR=PCR)
                
                if save_to_file:
                    if os.path.isfile(fn):
                        (count, inserts) = gac(fn, len(refseq),
                                                qual_min=qual_min,
                                                VERBOSE=VERBOSE)

                        counts_sample[PCR - 1] = count
                        count.dump(fn_out)

                        if VERBOSE >= 2:
                            print 'Allele counts saved:', samplename_pat, fragment

                else:
                    if os.path.isfile(fn_out):
                        count = np.load(fn_out)
                        counts_sample[PCR - 1] = count
                    elif os.path.isfile(fn):
                        (count, inserts) = gac(fn, len(refseq),
                                                qual_min=qual_min,
                                                VERBOSE=VERBOSE)
                        counts_sample[PCR - 1] = count

                

            counts.append(counts_sample)

