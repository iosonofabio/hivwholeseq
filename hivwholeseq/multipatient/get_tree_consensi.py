# vim: fdm=marker
'''
author:     Fabio Zanini
date:       22/06/14
content:    Take all consensi from sequenced samples connected to patient samples,
            make a multiple sequence alignment (MSA), and build a tree with it.

            Seq repetitions and PCR1/2 of the same sample should cluster, then within
            a patient. Else, contamination!
'''
# Modules
import os
import sys
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo

from hivwholeseq.patients.patients import load_patient, SamplePat
from hivwholeseq.patients.patients import load_samples_sequenced as lssp
from hivwholeseq.mapping_utils import align_muscle
from hivwholeseq.patients.filenames import get_consensi_alignment_filename, \
        get_consensi_tree_filename
from hivwholeseq.tree_utils import build_tree_fasttree
from hivwholeseq.reference import load_custom_reference


# Globals
refnames = ['38304', '38540', 'LAI-III']




# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Build tree of all patient samples')
    pats_or_samples = parser.add_mutually_exclusive_group(required=False)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to map (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--regions', nargs='+',
                        help='Regions to analyze (e.g. V3 F6)')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--no-summary', action='store_false', dest='summary',
                        help='Do not save results in a summary file')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    fragments = args.fragments
    submit = args.submit
    VERBOSE = args.verbose
    summary = args.summary


    refseq = load_custom_reference('HXB2', 'gb')

    # FIXME: finish to port this to JSON

    samples = lssp()
    # Collect all sequenced samples from patients
    if pnames is not None:
        samples = samples.loc[samples.patient.isin(pnames)]
        
    elif samplenames is not None:
        samples = samples.loc[samplenames]

    if VERBOSE >= 2:
        print 'samples', samples.index.tolist()

    for fragment in fragments:
        consensi = [load_custom_reference(refname+'_'+fragment) for refname in refnames]
        for samplename, sample in samples.iterrows():
            if VERBOSE >= 1:
                print samplename, fragment,
                if VERBOSE == 1:
                    print ''

            sample = SamplePat(sample)

            try:
                consensus = sample.get_consensus(fragment)
                consensus.id = sample.patient+'_'+sample.date.strftime('%Y-%m-%d')+'_'+samplename
                consensus.name = consensus.id
                consensi.append(consensus)
                if VERBOSE >= 2:
                    print 'OK'
            except IOError, OSError:
                if VERBOSE >= 2:
                    print 'ERROR: consensus file not found'
                continue

        if VERBOSE:
            print 'N consensi:', len(consensi)
                
        # Align
        if VERBOSE:
            print 'Align'
        ali = align_muscle(*consensi)
        ali_fn = get_consensi_alignment_filename('all', fragment)
        AlignIO.write(ali, ali_fn, 'fasta')

        # Make tree
        if VERBOSE:
            print 'Build tree'
        tree_fn = get_consensi_tree_filename('all', fragment)
        tree = build_tree_fasttree(ali_fn, VERBOSE=VERBOSE)
        Phylo.write(tree, tree_fn, 'newick')
