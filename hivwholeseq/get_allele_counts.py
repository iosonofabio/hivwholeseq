#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/08/13
content:    Get the allele frequencies out of a BAM file and a reference.
'''
# Modules
import os
import argparse
import subprocess as sp
import cPickle as pickle
from Bio import SeqIO

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.adapter_info import load_adapter_table
from hivwholeseq.filenames import get_mapped_filename, get_allele_counts_filename, \
        get_insert_counts_filename, get_coverage_filename, get_consensus_filename
from hivwholeseq.mapping_utils import convert_sam_to_bam
from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file
from hivwholeseq.fork_cluster import fork_get_allele_counts as fork_self
from hivwholeseq.samples import samples



# Globals
# Minimal quality required for a base to be considered trustful (i.e. added to 
# the allele counts), in phred score. Too high: lose coverage, too low: seq errors.
# Reasonable numbers are between 30 and 36.
qual_min = 35



# Functions
def get_allele_counts(data_folder, adaID, fragment, VERBOSE=0,
                      maxreads=1e10):
    '''Extract allele and insert counts from a bamfile'''

    # Read reference
    reffilename = get_consensus_filename(data_folder, adaID, fragment,
                                         trim_primers=True)
    refseq = SeqIO.read(reffilename, 'fasta')

    # Open BAM file
    # Note: the reads should already be filtered of unmapped stuff at this point
    bamfilename = get_mapped_filename(data_folder, adaID, fragment, type='bam',
                                      filtered=True)
    if not os.path.isfile(bamfilename):
        convert_sam_to_bam(bamfilename)

    # Call lower-level function
    return get_allele_counts_insertions_from_file(bamfilename, len(refseq),
                                                  qual_min=qual_min,
                                                  maxreads=maxreads, VERBOSE=VERBOSE)


def write_output_files(data_folder, adaID, fragment,
                       counts, inserts, coverage, VERBOSE=0):
    '''Write allele counts, inserts, and coverage to file'''
    if VERBOSE >= 1:
        print 'Write to file: '+adaID+' '+fragment

    # Save counts and coverage
    counts.dump(get_allele_counts_filename(data_folder, adaID, fragment))
    coverage.dump(get_coverage_filename(data_folder, adaID, fragment))

    # Convert inserts to normal nested dictionary for pickle
    inserts_dic = {k: dict(v) for (k, v) in inserts.iteritems()}
    with open(get_insert_counts_filename(data_folder, adaID, fragment), 'w') as f:
        pickle.dump(inserts_dic, f, protocol=-1)





# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Get allele counts')
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

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    submit = args.submit

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = dataset['adapters']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Iterate over all requested samples
    for adaID in adaIDs:

        # If the script is called with no fragment, iterate over all
        samplename = dataset['samples'][dataset['adapters'].index(adaID)]
        if not fragments:
            fragments_sample = [fr[:2] for fr in samples[samplename]['fragments']]
        else:
            from re import findall
            fragments_all = samples[samplename]['fragments']
            fragments_sample = []
            for fragment in fragments:
                frs = filter(lambda x: fragment in x, fragments_all)
                if len(frs):
                    fragments_sample.append(frs[0][:2])

        for fragment in fragments_sample:

            # Submit to the cluster self if requested
            if submit:
                fork_self(seq_run, adaID, fragment, VERBOSE=VERBOSE)
                continue

            # Get counts
            counts, inserts = get_allele_counts(data_folder, adaID, fragment,
                                                VERBOSE=VERBOSE)

            # Get coverage
            coverage = counts.sum(axis=1)

            # Save to file
            write_output_files(data_folder, adaID, fragment,
                               counts, inserts, coverage,
                               VERBOSE=VERBOSE)
