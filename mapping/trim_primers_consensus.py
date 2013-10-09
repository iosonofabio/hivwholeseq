# vim: fdm=marker
'''
author:     Fabio Zanini
date:       29/09/13
content:    Trim the inner PCR primers from a consensus sequence.

            Note: sometimes, the consensus sequences do not include the full
            primer!
'''
# Modules
import argparse
import numpy as np
from Bio import SeqIO
from Bio import pairwise2

from mapping.datasets import MiSeq_runs
from mapping.filenames import get_consensus_filename
from mapping.adapter_info import load_adapter_table
from mapping.primer_info import primers_inner



# Functions
def trim_primers(data_folder, adaID, fragment, subsample=False, VERBOSE=0):
    '''Trim the inner PCR primers from the consensus'''

    if VERBOSE >= 1:
        print 'Trimming primers from consensus:', adaID, fragment 

    # Resolve ambiguous fragment primers
    if fragment in ['F5a', 'F5b']:
        frag_gen = 'F5'
    else:
        frag_gen = fragment

    # Get the consensus
    reffilename = get_consensus_filename(data_folder, adaID, frag_gen,
                                         subsample=subsample) 
    refseq = SeqIO.read(reffilename, 'fasta')

    # Get the primers
    primer_fwd, primer_rev = primers_inner[fragment]

    # Use local alignment to trim the primers
    len_localali = 30

    # 1. FWD
    ref_fwd = str(refseq.seq)[:len_localali]
    ali_fwd = pairwise2.align.localms(ref_fwd, primer_fwd, 2, -1, -1.5, -0.1)[0][:2]
    # Check whether there is a chunk of primer at all
    alinp = np.array(map(list, ali_fwd))[:, -len(ali_fwd[0].lstrip('-')):len(ali_fwd[1].rstrip('-'))]
    frac_diff = (alinp[0] != alinp[1]).mean()
    if frac_diff > 0.3:
        # There is nothing to trim
        primer_fwd_end = 0
        if VERBOSE >= 1:
            'Fwd: nothing to trim'
    else:
        primer_fwd_end = len(ref_fwd) - ali_fwd[1][::-1].index(primer_fwd[-1])
        primer_fwd_end -= ali_fwd[0][:primer_fwd_end].count('-')

    # 2. REV
    ref_rev = str(refseq.seq)[-len_localali:]
    ali_rev = pairwise2.align.localms(ref_rev, primer_rev, 2, -1, -1.5, -0.1)[0][:2]
    # Check whether there is a chunk of primer at all
    alinp = np.array(map(list, ali_rev))[:, -len(ali_rev[1].lstrip('-')):len(ali_rev[0].rstrip('-'))]
    frac_diff = (alinp[0] != alinp[1]).mean()
    if frac_diff > 0.3:
        # There is nothing to trim
        primer_rev_start = len(refseq)
        if VERBOSE >= 1:
            'Rev: nothing to trim'
    else:
        primer_rev_start = len(refseq) - len(ref_rev) + ali_rev[1].index(primer_rev[0])
        primer_rev_start += ref_rev[len(ref_rev) - len(refseq) + primer_rev_start:].count('-')

    if VERBOSE >= 1:
        print 'Length of fwd primer in ref:', str(primer_fwd_end)+'/'+str(len(primer_fwd)), \
              'length of rev primer in ref:', str(len(refseq) - primer_rev_start)+'/'+str(len(primer_rev))

    # Trim
    refseq_new = refseq[primer_fwd_end: primer_rev_start]
    
    # Write output
    outputfilename = get_consensus_filename(data_folder, adaID, frag_gen,
                                         subsample=subsample, trim_primers=True)
    SeqIO.write(refseq_new, outputfilename, 'fasta')

    if VERBOSE:
        print 'Trimmed consensus written:', outputfilename



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Study minor allele frequency')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--subsample', action='store_true',
                        help='Apply only to a subsample of the reads')

    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    subsample = args.subsample

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = load_adapter_table(data_folder)['ID']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Iterate over all requested samples
    for fragment in fragments:
        for adaID in adaIDs:
            if fragment == 'F5':
                frag_orig = dataset['primerF5'][dataset['adapters'].index(adaID)]
            else:
                frag_orig = fragment
            trim_primers(data_folder, adaID, frag_orig, subsample=subsample,
                         VERBOSE=VERBOSE)
