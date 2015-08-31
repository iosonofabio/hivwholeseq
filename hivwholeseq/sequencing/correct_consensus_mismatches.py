# vim: fdm=marker
'''
author:     Fabio Zanini
date:       25/10/13
content:    Sometimes, an allele has a frequency very close to 0.5 and our initial
            subsample picks the wrong nucleotide. We correct these cases by using
            the full dataset.

            NOTE: in principle, we should do the same with insertions. That would,
            however, force us to remap, so for now ignore that.
'''
# Modules
import os
import shutil
import argparse
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO as SeqIO
import Bio.AlignIO as AlignIO

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.utils.miseq import alpha
from hivwholeseq.sequencing.primer_info import primers_inner as pri
from hivwholeseq.sequencing.adapter_info import load_adapter_table
from hivwholeseq.sequencing.filenames import get_consensus_filename, get_consensus_old_filename, \
        get_allele_counts_filename, get_coverage_filename
from hivwholeseq.sequencing.minor_allele_frequency import filter_nus




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
            cmat_af = alpha[nu.argmax(axis=0)]

            if len(cmat) != len(cmat_af):
                raise ValueError('The two consensi have a different length!')

            pos_diff = (cmat != cmat_af).nonzero()[0]

            # If they are the same, do nothing (we do not want useless backup files)
            if len(pos_diff) == 0:
                continue

            # Too many differences are a bad sign
            elif len(pos_diff) > 20:
                raise ValueError('The two consensi are too different')

            # Copy the consensus into a backup file
            output_old_filename =  get_consensus_old_filename(data_folder, adaID, fragment)
            output_filename =  get_consensus_filename(data_folder, adaID, fragment)
            if os.path.isfile(output_old_filename):
                os.remove(output_old_filename)
                if VERBOSE:
                    print 'Very old consensus removed'
            shutil.copy(output_filename, output_old_filename)
            if VERBOSE:
                print 'Consensus copied to back-up file'

            # Make the new consensus
            consensus_new = SeqRecord(Seq(cmat_af.tostring(), consensus.seq.alphabet),
                                      id=consensus.id, name=consensus.name,
                                      description=consensus.description)

            # Write to file
            SeqIO.write(consensus_new,
                        get_consensus_filename(data_folder, adaID, fragment),
                        'fasta')
            if VERBOSE:
                print 'New consensus written'
