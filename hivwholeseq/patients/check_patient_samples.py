# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/02/14
content:    Check how many and what samples we have sequenced from each patient.
'''
# Modules
import os
import argparse
from Bio import SeqIO
from Bio import AlignIO

from hivwholeseq.samples import samples as samples_seq
from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.filenames import get_merged_consensus_filename, get_consensus_filename
from hivwholeseq.patients.filenames import get_consensi_alignment_genomewide_filename, \
        get_consensi_alignment_filename
from hivwholeseq.patients.patients import patients
from hivwholeseq.mapping_utils import align_muscle



# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Check patient samples')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    VERBOSE = args.verbose

    # Get patient and the sequenced samples
    for p in patients:
        print p

