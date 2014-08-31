# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/03/14
content:    Load the initial reference for manual analysis.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO

from hivwholeseq.patients.patients import get_patient
from hivwholeseq.patients.filenames import get_initial_reference_filename
from hivwholeseq.sequence_utils import correct_genbank_features_load



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Load initial consensus')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--seqformat', default='fasta',
                        help='Specify a sequence format, e.g. gb')
    parser.add_argument('--save-genbank', dest='save_gb', action='store_true',
                        help='Save the annotated sequences in genbank format')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose
    seqformat = args.seqformat
    save_gb = args.save_gb

    patient = get_patient(pname)

    # If the script is called with no fragment, take genomewide reference
    if not fragments:
        fragments = ['genomewide']
    if VERBOSE >= 3:
        print 'fragments', fragments

    seqs = []
    for fragment in fragments:
        fn = get_initial_reference_filename(pname, fragment, format=seqformat)
        seq_rec = SeqIO.read(fn, seqformat)
        if seqformat in ['gb', 'genbank']:
            correct_genbank_features_load(seq_rec)
        seqs.append(seq_rec)

        if save_gb:
            if seqformat in ['gb', 'genbank']:
                print 'The input sequence is already in genbank format: skipping annotation and storage.'

            if VERBOSE:
                print 'Annotating sequence'
            from hivwholeseq.sequence_utils import annotate_sequence_genes
            annotate_sequence_genes(seq_rec, fragment=fragment, VERBOSE=VERBOSE)

            if VERBOSE:
                print 'Storing sequence in genbank format'
            from hivwholeseq.sequence_utils import correct_genbank_features_save
            correct_genbank_features_save(seq_rec)
            fngb = get_initial_reference_filename(pname, fragment, format='gb')
            SeqIO.write(seq_rec, fngb, 'gb')

