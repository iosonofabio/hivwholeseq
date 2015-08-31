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
from Bio import AlignIO

from hivwholeseq.patients.patients import load_patients, iterpatient
from hivwholeseq.utils.generic import write_json
from hivwholeseq.utils.tree import tree_to_json, correct_minimal_branches
from hivwholeseq.store.store_haplotypes_alignment_tree import extract_alignment
from hivwholeseq.website.filenames import get_haplotype_tree_filename
from hivwholeseq.website.filenames import get_haplotype_alignment_filename
from hivwholeseq.website.filenames import get_consensi_alignment_filename
from hivwholeseq.website import _regions_haplotypes as regions_all



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Store haplotypes and trees',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--regions', nargs='+', default=regions_all,
                        help='Regions to store (e.g. V3 PR)')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    for pname, patient in iterpatient(patients):
        for region in regions:
            print patient.code, patient.name, region

            print 'Get haplotype trees and alignments'
            tree = patient.get_local_tree(region)

            print 'Correct minimal-length branches to a smaller value'
            correct_minimal_branches(tree, VERBOSE=2)
            
            print 'Extract alignment from tree (so with duplicates)'
            ali_tree = extract_alignment(tree, VERBOSE=2)

            print 'Save alignment from tree, i.e. with duplicates'
            fn = get_haplotype_alignment_filename(patient.code, region, format='fasta')
            AlignIO.write(ali_tree, fn, 'fasta')

            # FIXME: select only leaves that are likely to exist (e.g. above
            # inverse template number etc.) for the tree, because we plot it!

            print 'Save tree (JSON)'
            fn = get_haplotype_tree_filename(patient.code, region, format='json')
            tree_json = tree_to_json(tree.root,
                                     fields=('DSI', 'sequence',
                                             'muts',
                                             'VL', 'CD4',
                                             'frequency',
                                             'count',
                                             'confidence'),
                                    )
            write_json(tree_json, fn, indent=1)


