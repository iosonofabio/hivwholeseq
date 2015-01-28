# vim: fdm=indent
'''
author:     Fabio Zanini
date:       06/08/14
content:    Check for contamination by building the matrix of minimal distance
            from any consensus for a bunch of reads.
'''
# Modules
import os
import sys
import argparse
from operator import itemgetter
import pysam
import numpy as np
from Bio import SeqIO
from matplotlib import cm
import matplotlib.pyplot as plt
from seqanpy import align_overlap

from hivwholeseq.sequencing.samples import load_sequencing_run, SampleSeq
from hivwholeseq.sequencing.filenames import get_consensus_filename, get_mapped_filename, \
        get_filter_mapped_summary_filename, get_mapped_suspicious_filename
from hivwholeseq.utils.mapping import get_ind_good_cigars, convert_sam_to_bam,\
        pair_generator, get_range_good_cigars



# Functions
def get_minimal_distance_hist(bamfilename, consensi, maxreads=1000, VERBOSE=0):
    '''Get histogram of minimal distance of reads from consensi'''

    conssi = map(''.join, consensi)
    m = np.zeros(len(consensi), int)
    n_good = 0
    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        for irp, reads in enumerate(pair_generator(bamfile)):
            if n_good == maxreads:
                break

            if VERBOSE >= 3:
                print n_good + 1, 'Checking mindist for:', reads[0].qname,

            # Assign names
            (read1, read2) = reads
            i_fwd = reads[0].is_reverse

            # Check a few things to make sure we are looking at paired reads
            if read1.qname != read2.qname:
                raise ValueError('Read pair '+str(irp)+': reads have different names!')

            # Ignore unmapped reads
            if read1.is_unmapped or read2.is_unmapped:
                if VERBOSE >= 2:
                    print 'Read pair '+read1.qname+': unmapped'
                continue
            
            # Ignore not properly paired reads (this includes mates sitting on
            # different fragments)
            if (not read1.is_proper_pair) or (not read2.is_proper_pair):
                if VERBOSE >= 2:
                    print 'Read pair '+read1.qname+': not properly paired'
                continue

            n_good += 1

            # Get all distances
            ds_pair = np.zeros_like(m)
            for ic, consensus in enumerate(consensi):
                conss = conssi[ic]
                dpair = 0
                for read in reads:
                    seq = read.seq
                    ali = align_overlap(conss, seq)
    
                    # NOTE: it is possible that we start before conss' start or end after
                    # its end, but that IS evidence that it's not contamination from there.
    
                    pos = conss.find(ali[1].replace('-', ''))
                    alim0 = np.fromstring(ali[1], 'S1')
                    alim1 = np.fromstring(ali[2], 'S1')
    
                    # Score subst
                    d = ((alim0 != alim1) & (alim0 != '-') & (alim1 != '-')).sum()
    
                    # Score insertions
                    gaps = alim0 == '-'
                    if gaps.sum():
                        n_gaps_borders = np.diff(gaps).sum()
                        n_gaps_borders += alim0[0] == '-'
                        n_gaps_borders += alim0[-1] == '-'
                        n_insertions = n_gaps_borders // 2
                        d += n_insertions
    
                    # Score deletions
                    gaps = alim1 == '-'
                    if gaps.sum():
                        n_gaps_borders = np.diff(gaps).sum()
                        n_gaps_borders -= alim1[0] == '-'
                        n_gaps_borders -= alim1[-1] == '-'
                        n_deletions = n_gaps_borders // 2
                        d += n_deletions
    
                    dpair += d
    
                ds_pair[ic] = dpair

                if VERBOSE >= 3:
                    print 'OK',

            m[ds_pair.argmin()] += 1
            if VERBOSE >= 3:
                print ''

    return m



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Check distance matrix from consensi',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. 28, 37)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help=('Verbosity level [0-3]'))
    parser.add_argument('--maxreads', type=int, default=1000,
                        help='Maximal number of reads to analyze')

    args = parser.parse_args()
    seq_run = args.run
    fragments = args.fragments
    VERBOSE = args.verbose
    maxreads = args.maxreads

    # Specify the dataset
    dataset = load_sequencing_run(seq_run)
    data_folder = dataset.folder

    samples = dataset.samples

    if not fragments:
        fragments = ['F'+str(i+1) for i in xrange(6)]

    matrices = {}
    for fragment in fragments:
        if VERBOSE:
            print fragment

        samples_frag = samples.loc[[os.path.isfile(SampleSeq(s).get_mapped_filename(fragment, type='bam', filtered=False))
                                   for sn, s in samples.iterrows()]]
        n_samples = len(samples_frag)

        consensi = [SeqIO.read(SampleSeq(s).get_consensus_filename(fragment), 'fasta') for sn, s in samples_frag.iterrows()]
        labels = [(sn, s.adapter) for sn, s in samples_frag.iterrows()]
        m = np.zeros((n_samples, n_samples), int)

        for si, (samplename, sample) in enumerate(samples_frag.iterrows()):
            if VERBOSE == 1:
                print samplename,
            elif VERBOSE >= 2:
                print samplename

            sample = SampleSeq(sample)
            adaID = sample.adapter
            bamfilename = get_mapped_filename(data_folder, adaID, fragment, type='bam',
                                              filtered=False)
            if not os.path.isfile(bamfilename):
                if VERBOSE >= 1:
                    print 'missing mapped file, skipping'
                continue

            m[si] = get_minimal_distance_hist(bamfilename, consensi, VERBOSE=VERBOSE, maxreads=maxreads)
    
        matrices[fragment] = m
