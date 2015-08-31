# vim: fdm=marker
'''
author:     Fabio Zanini
date:       20/11/13
content:    Check how many and what kind of reads crosses fragment boundaries.
'''
# Modules
import argparse
import numpy as np
from collections import Counter
from operator import itemgetter, attrgetter
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment as MSA
import pysam
from matplotlib import cm
import matplotlib.pyplot as plt

from mapping.miseq import alpha, alphal
from mapping.datasets import MiSeq_runs
from mapping.filenames import get_consensus_filename, get_mapped_filename, \
        get_allele_counts_filename, get_coverage_filename, get_divided_filenames
from mapping.mapping_utils import pair_generator



# Script
if __name__ == '__main__':

    # Parse input args: this is used to call itself recursively
    parser = argparse.ArgumentParser(description='Errors in HIV pure strains')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    fragments = args.fragments
    VERBOSE = args.verbose

    # Specify the dataset
    miseq_run = 28
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments


    ##########################################################################
    # Look for cross-fragment reads
    ##########################################################################
    # Focus on mix1
    adaID = 16

    ## We are going to have a special file for these, but for now filter the unmapped
    from mapping.primer_info import primers_coordinates_HXB2_inner as pcis
    from mapping.primer_info import primers_coordinates_HXB2_outer as pcos
    from mapping.primer_info import primers_inner, primers_outer
    from mapping.trim_and_divide import test_outer_primer
    from mapping.trim_and_divide import assign_to_fragment
    unmapped_filename = get_divided_filenames(data_folder, adaID, fragments)[-2]
    crossmapped_filename = get_divided_filenames(data_folder, adaID, fragments)[-3]
    fragments = ['F1', 'F2', 'F3', 'F4', 'F5b', 'F6']
    # This structure contains the fragment coordinates as of
    # - outer primers (row 0)
    # - inner primers (row 1)
    # - trimmed of all primers (row 2)
    frags_pos = np.zeros((3, 2, len(fragments)), int)
    for i, fragment in enumerate(fragments):
        pci = pcis[fragment]
        pco = pcos[fragment]
        frags_pos[0, :, i] = (pco[0][0], pco[1][1])
        frags_pos[1, :, i] = (pci[0][0], pci[1][1])
        frags_pos[2, :, i] = (pci[0][1], pci[1][0])
    # Since the reference is cropped, subtract from the positions F1 start
    # Note: now we are in the reference of the CROPPED HXB2, and start from 0!
    frags_pos -= frags_pos[1].min()

    # Make primers with masks for ambiguous nucleotides
    pr_outs = []
    for fragment in fragments:
        ptmps = primers_outer[fragment]
        for i, ptmp in enumerate(ptmps):
            ptmp = np.ma.array(list(ptmp), mask=[p not in alpha[:4] for p in ptmp])
            ptmps[i] = ptmp
        pr_outs.append(ptmps)

    with pysam.Samfile(unmapped_filename, 'rb') as input_file:
        with pysam.Samfile(crossmapped_filename, 'wb', template=input_file) as output_file:
            for reads in pair_generator(input_file):
                # Truly unmapped stuff we do not want
                if reads[0].is_unmapped or (not reads[0].is_proper_pair):
                    if VERBOSE >= 3:
                        print 'Read pair unmapped/unpaired/tiny:', reads[0].qname
                    continue

                # Stuff from the outer primers we do not want
                if test_outer_primer(reads, pr_outs, frags_pos[1, 1, -1]):
                    if VERBOSE >= 3:
                        print 'Read pair from outer primer:', reads[0].qname
                    continue

                # Check fragment assignment
                frags_pair = assign_to_fragment(reads, frags_pos[1])
                if len(frags_pair) > 0:
                    continue

                # The other reads are good
                output_file.write(reads[0])
                output_file.write(reads[1])


