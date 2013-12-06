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



# Globals
# Cluster submit
import hivwholeseq
JOBDIR = hivwholeseq.__path__[0].rstrip('/')+'/'
JOBSCRIPT = JOBDIR+'filter_allele_frequencies.py'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
cluster_time = '0:59:59'
vmem = '2G'





# Functions
def fork_self(miseq_run, adaID, fragment, VERBOSE=3):
    '''Fork self for each adapter ID and fragment'''
    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'acn '+'{:02d}'.format(adaID)+' '+fragment,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', miseq_run,
                 '--adaIDs', adaID,
                 '--fragments', fragment,
                 '--verbose', VERBOSE,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def filter_nus(counts, coverage):
    '''Filter allele frequencies from the four read types'''
    pvals = np.zeros((len(alpha), counts.shape[-1]))

    # Divide binarily
    nocounts = (coverage - counts.swapaxes(0, 1)).swapaxes(0, 1)

    # Sum read1 and read2
    counts_f = counts[0] + counts[2]
    counts_b = counts[1] + counts[3]
    nocounts_f = nocounts[0] + nocounts[2]
    nocounts_b = nocounts[1] + nocounts[3]
    nus_f = 1.0 * counts_f / (coverage[0] + coverage[2] + 1)
    nus_b = 1.0 * counts_b / (coverage[1] + coverage[3] + 1)

    # Test chi^2 for consistency across read 1 and read2,
    # using pseudocounts
    from scipy.stats import chi2_contingency
    for j in xrange(len(alpha)):
        for i in xrange(counts.shape[-1]):
            cm = np.array([[counts_f[j, i], nocounts_f[j, i]],
                           [counts_b[j, i], nocounts_b[j, i]]], int)
            chi2, pval = chi2_contingency(cm + 1)[:2]
            pvals[j, i] = pval

    errs = zip(*(((pvals < 1e-6) & (np.abs(nus_f - nus_b) > 1e-4)).nonzero()))
    errs.sort(key=itemgetter(1))

    # Take the mean of the two fed and rev for non-errors
    nu_filtered = np.ma.masked_all((len(alpha), counts.shape[-1]))
    covtot = coverage.sum(axis=0)
    ind = covtot > 0
    # Geometric mean? Not much changes
    #nu_filtered[:, ind] = 1.0 * (counts_f + counts_b)[:, ind] / covtot[ind]
    nu_filtered[:, ind] = np.sqrt(nus_f * nus_b)[:, ind]

    # Insert corrections
    if VERBOSE:
        print '{:4s}'.format('pos'), '{:3s}'.format('nuc'), \
                '{:5s}'.format('cov fw'), '{:10s}'.format('nu fw'), \
                '{:10s}'.format('nu re'), '{:5s}'.format('cov rv')
    for ai, pos in errs:
        # Take the allele farthest away from 0.5
        # (the sum to 1 is not guaranteed!)
        nu_tmp = np.array([nus_f[ai, pos], nus_b[ai, pos]])
        ind_tmp = np.argmax(np.abs(nu_tmp - 0.5))
        nu_filtered[ai, pos] = nu_tmp[ind_tmp]

        if VERBOSE and ((nus_f[ai, pos] < 0.5) or (nus_b[ai, pos] < 0.5)):
            print '{:4d}'.format(pos), alpha[ai], \
                    '{:7d}'.format(coverage[[0, 2], pos].sum()), \
                    '{:1.1e}'.format(nus_f[ai, pos]), \
                    '{:1.1e}'.format(nus_b[ai, pos]), \
                    '{:7d}'.format(coverage[[1, 3], pos].sum())

    # Get rid of the mask if not needed
    if not nu_filtered.mask.any():
        nu_filtered = nu_filtered.data

    return nu_filtered


def write_frequencies(data_folder, adaID, fragment, nu_filtered, VERBOSE=0):
    '''Write the corrected allele frequencies to file'''
    if VERBOSE:
        print 'Storing allele frequencies to file:', adaID, fragment
    nu_filtered.dump(get_allele_frequencies_filename(data_folder, adaID, fragment))



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
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')

    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    submit = args.submit

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
    for adaID in adaIDs:
        for fragment in fragments:

            # Submit to the cluster self if requested
            if submit:
                fork_self(data_folder, adaID, fragment, VERBOSE=VERBOSE)
                continue

            # Get coverage and counts
            counts = np.load(get_allele_counts_filename(data_folder, adaID, fragment))
            coverage = np.load(get_coverage_filename(data_folder, adaID, fragment))
    
            # Filter the minor frequencies by comparing the read types
            nu_filtered = filter_nus(counts, coverage)

            # Write output
            write_frequencies(data_folder, adaID, fragment, nu_filtered,
                              VERBOSE=VERBOSE)

