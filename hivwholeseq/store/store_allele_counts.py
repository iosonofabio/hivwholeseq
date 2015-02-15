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

from hivwholeseq.patients.samples import load_samples_sequenced as lssp
from hivwholeseq.patients.samples import SamplePat
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_mapped_filtered_filename, get_allele_counts_filename
from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file as gac
from hivwholeseq.cluster.fork_cluster import fork_get_allele_counts_patient as fork_self 



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Get allele counts',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=True)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to analyze (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--fragments', nargs='+',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save the allele counts to file')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--plot', action='store_true',
                        help='Plot the allele frequency trajectories')
    parser.add_argument('--logit', action='store_true',
                        help='use logit scale (log(x/(1-x)) in the plots')
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
    use_plot = args.plot

    if use_plot:
        import matplotlib.pyplot as plt

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

    counts_all = []
    for fragment in fragments:
        counts = []
        for samplename, sample in samples.iterrows():
            if submit:
                fork_self(samplename, fragment, VERBOSE=VERBOSE, qual_min=qual_min)
                continue

            if VERBOSE >= 1:
                print fragment, samplename

            sample = SamplePat(sample)
            pname = sample.patient
            refseq = SeqIO.read(get_initial_reference_filename(pname, fragment), 'fasta')

            fn_out = sample.get_allele_counts_filename(fragment, PCR=PCR, qual_min=qual_min)
            fn = sample.get_mapped_filtered_filename(fragment, PCR=PCR, decontaminated=True) #FIXME
            
            # If --save, recalculate and save
            if save_to_file:
                if os.path.isfile(fn):
                    (count, inserts) = gac(fn, len(refseq),
                                            qual_min=qual_min,
                                            VERBOSE=VERBOSE)

                    count.dump(fn_out)

                    if VERBOSE >= 2:
                        print 'Allele counts saved:', samplename, fragment

                    counts.append(count)

            # If not --save, try to load from file and recalculate as a fallback
            elif os.path.isfile(fn_out):
                count = np.load(fn_out)
                counts.append(count)

            elif os.path.isfile(fn):
                (count, inserts) = gac(fn, len(refseq),
                                        qual_min=qual_min,
                                        VERBOSE=VERBOSE)
                counts.append(count)


            if use_plot:
                cou = count.sum(axis=0)
                x = np.tile(np.arange(cou.shape[1]), (cou.shape[0], 1))
                color = np.tile(np.arange(cou.shape[0]), (cou.shape[1], 1)).T

                fig, ax = plt.subplots(figsize=(12, 6))
                
                ax.scatter(x, cou + 0.1, lw=2, c=color)
                ax.set_xlabel('Position [bp]')
                ax.set_ylabel('Coverage')
                ax.set_xlim(-1, cou.shape[-1])
                ax.set_ylim(ymin=0.09)
                ax.set_yscale('log')
                ax.grid(True)
                ax.set_title(samplename)

    if use_plot:
        plt.ion()
        plt.show()
