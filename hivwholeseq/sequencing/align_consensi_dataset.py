# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/02/14
content:    Align the fragment and, if present, genomewide consensi in the dataset
            for qualitative assesment.
'''
# Modules
import os
import argparse
from collections import defaultdict
from Bio import SeqIO
import numpy as np

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.sequencing.filenames import get_consensus_filename, get_merged_consensus_filename, \
        get_consensi_alignment_dataset_filename
from hivwholeseq.sequencing.samples import samples
from hivwholeseq.utils.mapping import align_muscle



# Functions
def align_consensi_dataset(dataset, adaIDs, fragments, VERBOSE=0):
    '''Align consensi from different samples in a dataset'''

    data_folder = dataset['folder']

    # Collect consensi
    if VERBOSE >= 1:
        print 'Collecting consensi...',
    consensi = defaultdict(dict)
    for adaID in adaIDs:
        samplename = dataset['samples'][dataset['adapters'].index(adaID)]
        fragments_sample = samples[samplename]['fragments']
        for frag in fragments_sample:
            frag_gen = frag[:2]
            if frag_gen not in fragments:
                continue
            con_fn = get_consensus_filename(data_folder, adaID, frag_gen)
            if os.path.isfile(con_fn):
                con = SeqIO.read(con_fn, 'fasta')
                consensi[frag_gen][adaID] = con

        if 'genomewide' in fragments:
            frag_gens = [frag[:2] for frag in fragments_sample]
            con_gw_fn = get_merged_consensus_filename(data_folder, adaID, frag_gens)
            if os.path.isfile(con_gw_fn):
                con = SeqIO.read(con_gw_fn, 'fasta')
                consensi['genomewide'][adaID] = con
    
    if VERBOSE >= 1:
        print 'done.'
        print 'Aligning...',

    # Align
    alis = {}
    for (frag, con_dict) in consensi.iteritems():
        if VERBOSE >= 2:
            print frag,
        ali_frag = align_muscle(*(con_dict.values()))
        alis[frag] = ali_frag

    if VERBOSE >= 1:
        print 'done.'

    return alis



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Align consensi from a dataset')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = dataset['adapters']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    if fragments is None:
        fragments = ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'genomewide']

    # Iterate over all requested samples
    alis = align_consensi_dataset(dataset, adaIDs, fragments, VERBOSE=VERBOSE)
    for (frag, ali) in alis.iteritems():
        SeqIO.write(ali, get_consensi_alignment_dataset_filename(data_folder, frag),
                    'fasta')
