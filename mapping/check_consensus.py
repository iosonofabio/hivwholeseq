# vim: fdm=marker
'''
author:     Fabio Zanini
date:       24/09/13
content:    Check iterative consensus with reference (if appropriate) and between
            fragments at overlaps.
'''
# Modules
import os
import sys
import subprocess as sp
import time
import argparse
import re
from operator import *
from itertools import izip
from collections import defaultdict
from collections import Counter
import numpy as np
import pysam
import Bio.SeqIO as SeqIO
import Bio.AlignIO as AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2

# Horizontal import of modules from this folder
from mapping.adapter_info import load_adapter_table, foldername_adapter
from mapping.miseq import alpha, read_types
from mapping.mapping_utils import stampy_bin, subsrate, convert_sam_to_bam,\
        pair_generator, get_ind_good_cigars
from mapping.filenames import get_HXB2_fragmented, get_read_filenames,\
        get_HXB2_index_file, get_HXB2_hash_file, get_consensus_filename



# Globals
VERBOSE = 3
# FIXME
from mapping.datasets import dataset_2 as dataset
data_folder = dataset['folder']
references = {2: 'NL4-3', 4: 'SHIV_SF162'}


# Functions
def find_fragment(refseq, seq):
    '''Find a fragment location in a sequece (do not align)'''
    refs = str(refseq.seq)
    s = str(seq.seq)

    # Align the seq back to the reference
    for bsize in [40, 30, 20, 15, 10, 7]:
        i = 0
        while i + bsize < len(s):
            start = refs.find(s[i:i+bsize])
            if start != -1:
                start -= i
                break
            i += bsize
        if start != -1:
            break
    if start == -1:
        return None

    for bsize in [40, 30, 20, 15, 10, 7]:
        i = len(seq)
        while i - bsize > 0:
            end = refs[start:].find(s[i - bsize: i])
            if end != -1:
                end += bsize + len(seq) - i + start
                break
            i -= bsize
        if end != -1:
            break

    refsf = refseq[start:end]
    return refsf


def find_overlap(seq1, seq2):
    '''Find the overlap between two consecutive fragments'''
    s1 = str(seq1.seq)
    s2 = str(seq2.seq)

    #FIXME: invert the bsize and while loops!

    # Align the seq back to the reference
    for bsize in [40, 30, 20, 15, 10, 7]:
        i = 0
        while i + bsize < len(s2):
            start = s1.find(s2[i:i+bsize])
            if start != -1:
                start -= i
                break
            i += bsize
        if start != -1:
            break

    for bsize in [40, 30, 20, 15, 10, 7]:
        i = len(s1)
        while i - bsize > 0:
            end = s2.find(s1[i - bsize: i])
            if end != -1:
                end += bsize + len(s1) - i 
                break
            i -= bsize
        if end != -1:
            break

    # Check whether anything survived
    if (start == -1) or (end == -1):
        return None

    return (seq1[start:], seq2[:end])


def align_via_muscle(refseq, seq):
    '''Align two sequences via MUSCLE'''
    from Bio.Align.Applications import MuscleCommandline
    muscle_cline = MuscleCommandline(diags=True)
    import subprocess as sp
    import sys
    child = sp.Popen(str(muscle_cline),
                     stdin=sp.PIPE,
                     stdout=sp.PIPE,
                     stderr=sp.PIPE,
                     shell=True)
    SeqIO.write([refseq, seq], child.stdin, "fasta")
    child.stdin.close()
    child.stderr.close()
    align = AlignIO.read(child.stdout, "fasta")
    child.stdout.close()
    return align



# Script
if __name__ == '__main__':

    # Parse input args: this is used to call itself recursively
    parser = argparse.ArgumentParser(description='Map HIV reads recursively')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=5,
                        help='Verbosity level [0-3]')
    parser.add_argument('--reference', default='',
                        help='Also compare to a reference (e.g. NL4-3)')

    args = parser.parse_args()
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    reference = args.reference

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

    for adaID in adaIDs:

        ## If there is a reference, check it
        #if (adaID in [2, 4]) and (not reference):
        #    reference = references[adaID]
        if reference:
            refseq = SeqIO.read(data_folder+'reference/'+reference+'.fasta', 'fasta')
            
            for fragment in fragments:
                seq = SeqIO.read(get_consensus_filename(data_folder, adaID, fragment),
                                 'fasta')
                refsf = find_fragment(refseq, seq)

                # Align
                align = align_via_muscle(refsf, seq)

                # Look for mutations
                ali = np.array(align)
                muts = (ali[1] != ali[0]).nonzero()[0]
                if len(muts):
                    print adaID, fragment, muts

                #FIXME
                #import ipdb; ipdb.set_trace()

        # Look at the overlapping regions
        for fragment1 in fragments:
            fragment2 = 'F'+str(int(fragment1[-1])+1)
            if fragment2 not in fragments:
                continue

            seq1 = SeqIO.read(get_consensus_filename(data_folder, adaID, fragment1),
                              'fasta')
            seq2 = SeqIO.read(get_consensus_filename(data_folder, adaID, fragment2),
                              'fasta')

            # Find overlap
            (seqo1, seqo2) = find_overlap(seq1, seq2)

            # Align
            align = align_via_muscle(seqo1, seqo2)

            # Look for mutations
            ali = np.array(align)
            muts = (ali[1] != ali[0]).nonzero()[0]
            if len(muts):
                print adaID, fragment1+' - '+fragment2, ali.shape[1], muts








