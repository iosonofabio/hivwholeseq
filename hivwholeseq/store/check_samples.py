# vim: fdm=marker
'''
author:     Fabio Zanini
date:       15/12/14
content:    Check status of single samples in the pipeline.
'''
# Modules
import os
import argparse
import datetime
from hivwholeseq.patients.patients import load_samples_sequenced, Patient, \
        SamplePat
from hivwholeseq.patients.filenames import get_mapped_filtered_filename
from hivwholeseq.utils.generic import modification_date
from hivwholeseq.store.check_patients import (
    pretty_print_info, pretty_print_info_genomewide)


# Globals
title_len = 15
cell_len = 7
get_decontaminated_filename = lambda *args, **kwargs: get_mapped_filtered_filename(*args, decontaminated=True, **kwargs)



# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Check patient samples')
    parser.add_argument('--samples', nargs='+',
                        help='Samples to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    VERBOSE = args.verbose
    samplenames = args.samples

    samples = load_samples_sequenced()
    if samplenames is not None:
        samples = samples.loc[samplenames]

    for samplename, sample in samples.iterrows():
        sample = SamplePat(sample)

        mod_dates = {}

        pretty_print_info(sample, 'Map + filter', 'filter',
                          'get_mapped_filtered_filename',
                          None,#'reference',
                          mod_dates,
                          VERBOSE=VERBOSE)

        pretty_print_info(sample, 'Decontaminate', 'decontaminate',
                          get_decontaminated_filename,
                          'filter', mod_dates,
                          VERBOSE=VERBOSE)

        pretty_print_info(sample, 'Consensus', 'consensus',
                          'get_consensus_filename',
                          'decontaminate', mod_dates,
                          VERBOSE=VERBOSE)

        pretty_print_info_genomewide(sample, 'Cons genomewide', 'consensus',
                                     'get_consensus_filename',
                                     mod_dates,
                                     VERBOSE=VERBOSE)

        pretty_print_info(sample, 'Allele counts', 'allele counts',
                          'get_allele_counts_filename',
                          'decontaminate', mod_dates,
                          VERBOSE=VERBOSE)

        pretty_print_info(sample, 'Allele cocounts', 'allele cocounts',
                          'get_allele_cocounts_filename',
                          'decontaminate', mod_dates,
                          VERBOSE=VERBOSE)


