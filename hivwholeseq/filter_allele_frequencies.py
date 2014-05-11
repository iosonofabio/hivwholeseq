#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/08/13
content:    Correct the allele frequencies comparing read types and write to file.
'''
# Modules
import subprocess as sp
import argparse
from operator import itemgetter
import numpy as np

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.miseq import alpha
from hivwholeseq.filenames import get_allele_counts_filename, get_coverage_filename, \
        get_allele_frequencies_filename
from hivwholeseq.adapter_info import load_adapter_table
from hivwholeseq.one_site_statistics import filter_nus, plot_SFS_folded
from hivwholeseq.fork_cluster import fork_filter_allele_frequencies as fork_self
from hivwholeseq.samples import samples



# Functions
def write_frequency_files(data_folder, adaID, fragment, nu_filtered, VERBOSE=0):
    '''Write the corrected allele frequencies to file'''
    if VERBOSE:
        print 'Storing allele frequencies to file:', adaID, fragment

    nu_filtered.dump(get_allele_frequencies_filename(data_folder, adaID, fragment))



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Study minor allele frequency')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--no-summary', action='store_false', dest='summary',
                        help='Do not save results in a summary file')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    submit = args.submit
    summary = args.summary

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = load_adapter_table(data_folder)['ID']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # Iterate over all requested samples
    for adaID in adaIDs:

        # If the script is called with no fragment, iterate over all
        samplename = dataset['samples'][dataset['adapters'].index(adaID)]
        if fragments is None:
            fragments_sample = [fr[:2] for fr in samples[samplename]['fragments']]
        else:
            fragments_sample = fragments

        if VERBOSE >= 3:
            print 'adaID:', adaID+', fragments:', fragments_sample

        for fragment in fragments_sample:

            # Submit to the cluster self if requested
            if submit:
                fork_self(data_folder, adaID, fragment, VERBOSE=VERBOSE,
                          summary=summary)
                continue

            # Get coverage and counts
            counts = np.load(get_allele_counts_filename(data_folder, adaID, fragment))
            if len(counts.shape) == 2:
                import warnings
                warnings.warn('Counts not divided by read type: will normalize instead of filter!')
                nu_filtered = 1.0 * counts / counts.sum(axis=0)
    
            else:
                # Filter the minor frequencies by comparing the read types
                nu_filtered = filter_nus(counts)

            # Write output
            write_frequency_files(data_folder, adaID, fragment, nu_filtered,
                                  VERBOSE=VERBOSE)

            if summary:
                import matplotlib.pyplot as plt
                was_interactive = plt.isinteractive()
                plt.ioff()
                plot_SFS_folded(data_folder, adaID, fragment, nu_filtered,
                                VERBOSE=VERBOSE, savefig=True)
                plt.interactive(was_interactive)
