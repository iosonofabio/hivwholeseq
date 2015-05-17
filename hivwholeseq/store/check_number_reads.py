# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/05/15
content:    Check the number of reads for samples and patients.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO
import pandas as pd

from hivwholeseq.patients.samples import load_samples_sequenced as lssp
from hivwholeseq.patients.samples import itersample
from hivwholeseq.utils.argparse import PatientsAction


# Functions
def collect_nreads(sample, fragments, VERBOSE=0, title_len=15, cell_len=9):
    '''Collect number of reads from each sample'''
    n_reads = []

    title = sample.name
    line = ('{:<'+str(title_len)+'}').format(title+':')
    for fragment in fragments:
        try:
            n = sample.get_number_reads(fragment)
        except IOError:
            n = 0

        if n > 1000:
            nfmt = str(n/1000)+'k'
        else:
            nfmt = str(n)

        line = line+fragment+': '+\
            ('{:>'+str(cell_len - len(fragment) - 1)+'}').format(nfmt)+'  '

        n_reads.append({'pname': sample.patient,
                        'samplename': sample.name,
                        'fragment': fragment,
                        'n': n})

    if VERBOSE >= 1:
        print line

    return pd.DataFrame(n_reads)



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Get allele counts',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=False)
    pats_or_samples.add_argument('--patients', action=PatientsAction,
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to analyze (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--fragments', nargs='+',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    fragments = args.fragments
    VERBOSE = args.verbose

    samples = lssp()
    if pnames is not None:
        samples = samples.loc[samples.patient.isin(pnames)]
    elif samplenames is not None:
        samples = samples.loc[samples.index.isin(samplenames)]

    if VERBOSE >= 3:
        print 'samples', samples.index.tolist()

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    n_reads = []
    for samplename, sample in itersample(samples):
        n_reads.append(collect_nreads(sample, fragments, VERBOSE=VERBOSE))
    n_reads = pd.concat(n_reads)
