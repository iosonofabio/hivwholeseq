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

from hivwholeseq.utils.argparse import PatientsAction
from hivwholeseq.patients.samples import load_samples_sequenced as lssp
from hivwholeseq.patients.samples import SamplePat
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_mapped_filtered_filename, get_allele_counts_filename
from hivwholeseq.one_site_statistics import get_allele_counts_aa_from_file as gac
from hivwholeseq.cluster.fork_cluster import fork_get_allele_counts_aa_patient as fork_self 



# Functions




# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get allele counts of amino acids',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=False)
    pats_or_samples.add_argument('--patients', action=PatientsAction,
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to analyze')
    parser.add_argument('--proteins', nargs='+', required=True,
                        help='Proteins to analyze (e.g. PR IN)')
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
    proteins = args.proteins
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

    counts_all = []
    for protein in proteins:
        counts = []
        for samplename, sample in samples.iterrows():
            if submit:
                fork_self(samplename, protein, VERBOSE=VERBOSE, qual_min=qual_min)
                continue

            if VERBOSE >= 1:
                print protein, samplename

            sample = SamplePat(sample)

            # NOTE: How do we find what fragment covers the protein? Well, a
            # protein can happily cross fragments. Since each
            # codon is independent, we should iterate over codons. We do not
            # do that for efficiency reasons. Instead, we identify all potential
            # fragments and split the protein into full codon chunks covered by
            # a single fragment.
            fragment_rois = sample.get_fragments_covered(protein,
                                                         include_coordinates=True)
            
            refseq = sample.get_reference(protein)
            fn_out = sample.get_allele_counts_filename(protein, PCR=PCR,
                                                       qual_min=qual_min,
                                                       type='aa')

            from hivwholeseq.utils.sequence import alphaa
            count = np.zeros((len(alphaa), len(refseq) // 3), int)
            for frroi in fragment_rois:
                fragment = frroi['name']
                start_fr, end_fr = frroi['fragment']
                start, end = frroi['roi']

                # Check that we align with codons
                rf = start % 3
                if rf:
                    start_fr += 3 - rf
                    start += 3 - rf
                rf = end % 3
                if rf:
                    end_fr -= rf
                    end -= rf
                               
                fn = sample.get_mapped_filtered_filename(fragment, PCR=PCR,
                                                         decontaminated=True)
                
                if not os.path.isfile(fn):
                    if VERBOSE >= 2:
                        print 'SKIP'
                    continue
                
                if VERBOSE >= 2:
                    print 'Get allele counts for amino acids'
                cou = gac(fn, start_fr, end_fr, qual_min=qual_min,
                          VERBOSE=VERBOSE)
                # We do not care about fwd/reverse
                count[:, start // 3: end // 3] = cou.sum(axis=0)

            counts.append(count)

            if save_to_file:
                if VERBOSE >= 2:
                    print 'Save allele counts:', samplename, protein
                count.dump(fn_out)

            if use_plot:
                if VERBOSE >= 2:
                    print 'Plot'
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
