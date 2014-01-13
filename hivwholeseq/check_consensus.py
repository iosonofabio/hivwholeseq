# vim: fdm=marker
'''
author:     Fabio Zanini
date:       24/09/13
content:    Check whether the subsample consensus is also the real consensus of
            the full mapped dataset.
'''
# Modules
import argparse
import numpy as np
import Bio.SeqIO as SeqIO
import Bio.AlignIO as AlignIO

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.miseq import alpha
from hivwholeseq.adapter_info import load_adapter_table
from hivwholeseq.filenames import get_consensus_filename, \
        get_allele_counts_filename, get_coverage_filename
from hivwholeseq.minor_allele_frequency import filter_nus



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Check consensus')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--reference', default='',
                        help='Also compare to a reference (e.g. NL4-3)')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    reference = args.reference

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = MiSeq_runs[seq_run]['adapters']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Iterate over samples and fragments
    for adaID in adaIDs:
        for fragment in fragments:
            consensus = SeqIO.read(get_consensus_filename(data_folder, adaID, fragment),
                                   'fasta')
            cmat = np.array(consensus)

            counts = np.load(get_allele_counts_filename(data_folder, adaID, fragment))
            coverage = np.load(get_coverage_filename(data_folder, adaID, fragment))
            nu = filter_nus(counts, coverage, VERBOSE=VERBOSE)

            # Note: not-covered positions are filtered, but argmax cannot work
            # with masked arrays
            cmat_af = alpha[nu.argmax(axis=0)]
            if hasattr(nu, 'mask'):
                cmat_af[nu.mask.all(axis=0)] = 'N'

            # Check for consistency first
            if len(cmat) != len(cmat_af):
                print 'Consensus has a different length from allele frequency \
                        matrix... WTF?'

            # Do not actually align, it makes a huge mess (we miss mistakes)
            ali = [cmat, cmat_af]
            if (ali[0] != ali[1]).any():
                print 'adaID', adaID, fragment, 'consensus not self-coherent:'
                print 'differences (cons -> cons_af),'
                for pos in (ali[0] != ali[1]).nonzero()[0]:
                    print pos, ali[0][pos], '->', ali[1][pos], \
                            ali[1][max(0, pos-10): pos +1].tostring(), \
                            counts[:, (alpha == ali[0][pos]).nonzero()[0][0], pos], \
                            counts[:, (alpha == ali[1][pos]).nonzero()[0][0], pos]
