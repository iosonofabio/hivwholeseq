# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/04/14
content:    Try to guess recombination from the reference sample of subtype C,
            using the internal diversity.
'''
# Modules
import os
import sys
import argparse
from operator import itemgetter
from collections import defaultdict
import numpy as np
import pysam
from Bio import SeqIO

from hivwholeseq.miseq import alpha
from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.sequencing.filenames import get_consensus_filename, \
        get_allele_frequencies_filename, \
        get_mapped_filename
from hivwholeseq.mapping_utils import align_muscle, pair_generator


# Globals
sd = {'run': 'Tue59', 'adaID': 'N5-S4'}



# Script
if __name__ == '__main__':


    # Parse input args: this is used to call itself recursively
    parser = argparse.ArgumentParser(description='Errors in HIV pure strains')
    parser.add_argument('--fragment', default='F2',
                        help='Fragment to analyze (e.g. F2)')
    parser.add_argument('--maxreads', type=int, default=100,
                        help='Number of reads analyzed')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    fragment = args.fragment
    maxreads = args.maxreads
    VERBOSE = args.verbose


    dataset = MiSeq_runs[sd['run']]
    data_folder = dataset['folder']
    adaID = sd['adaID']

    cons = SeqIO.read(get_consensus_filename(data_folder, adaID, fragment), 'fasta')
    consm = np.array(cons)
    af = np.load(get_allele_frequencies_filename(data_folder, adaID, fragment))

    # Get minor alleles
    af2 = np.zeros(af.shape[-1])
    a2 = np.zeros(af.shape[-1], 'S1')
    ai2 = np.zeros(af.shape[-1], int)
    for i, afi in enumerate(af.T):
        j = np.argsort(afi)[-2]
        af2[i] = afi[j]
        a2[i] = alpha[j]
        ai2 = j

    # Get minor alleles at frequency ~0.5
    pos_poly = (af2 > 0.3).nonzero()[0]
    print 'Popular polymorphic positions:', pos_poly

    # Pick a group of those
    pos_poly = pos_poly[5: 9]
    a_poly = a2[pos_poly]
    a_anc = consm[pos_poly]
    start = pos_poly[0]
    end = pos_poly[-1] + 1


    # Go down to the reads and wonder about linkage
    hc = np.zeros([2 for p in pos_poly], int)
    n_cov = 0
    bamfilename = get_mapped_filename(data_folder, adaID, fragment, filtered=True)
    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        for irp, read_pair in enumerate(pair_generator(bamfile)):
            if (VERBOSE >= 2) and not ((irp + 1) % 1000):
                print irp + 1

            if irp == maxreads:
                if VERBOSE >= 2:
                    print 'Maxreads reached:', maxreads
                break

            i_fwd = read_pair[0].is_reverse
            i_rev = not i_fwd

            start_read_fwd = read_pair[i_fwd].pos
            end_read_fwd = start_read_fwd + \
                    sum(bl for (bt, bl) in read_pair[i_fwd].cigar if bt in (0, 2))
            start_read_rev = read_pair[i_rev].pos
            end_read_rev = start_read_rev + \
                    sum(bl for (bt, bl) in read_pair[i_rev].cigar if bt in (0, 2))

            covered = np.zeros(len(pos_poly), bool)
            covered[(pos_poly > start_read_fwd) & (pos_poly < end_read_fwd)] = True
            covered[(pos_poly > start_read_rev) & (pos_poly < end_read_rev)] = True

            if (not covered.all()) or (start_read_rev > end_read_fwd):
                continue

            markers_pair = []
            for read in read_pair:
                markers = {}
                pos_ref = read.pos
                pos_read = 0
                for (bt, bl) in read.cigar:
                    if bt == 1:
                        pos_read += bl
                    elif bt == 2:
                        ind_mark_pot = ((pos_poly >= pos_ref) & (pos_poly < pos_ref + bl)).nonzero()[0]
                        for i in ind_mark_pot:
                            if a_anc[i] == '-':
                                markers[pos_poly[i]] = 'A'
                            elif a_poly[i] == '-':
                                markers[pos_poly[i]] = 'M'
                            else:
                                markers[pos_poly[i]] ='X'

                        pos_ref += bl

                    else:
                        ind_mark_pot = ((pos_poly >= pos_ref) & (pos_poly < pos_ref + bl)).nonzero()[0]
                        for i in ind_mark_pot:
                            nuc = read.seq[pos_read + pos_poly[i] - pos_ref]
                            if nuc == a_anc[i]:
                                markers[pos_poly[i]] = 'A'
                            elif nuc == a_poly[i]:
                                markers[pos_poly[i]] = 'M'
                            else:
                                markers[pos_poly[i]] ='X'

                        pos_ref += bl
                        pos_read += bl

                markers_pair.append(markers)

            # Consensus within the read: exclude markers that are not agreed upon
            markers_pair_new = {}
            for pos in markers_pair[0]:
                if (pos not in markers_pair[1]) or (markers_pair[0][pos] == markers_pair[1][pos]):
                    markers_pair_new[pos] = markers_pair[0][pos]
            for pos in markers_pair[1]:
                if (pos not in markers_pair[0]):
                    markers_pair_new[pos] = markers_pair[1][pos]
            markers_pair_new = list(markers_pair_new.iteritems())
            markers_pair_new.sort(key=itemgetter(0))

            if len(markers_pair_new) < len(pos_poly):
                continue

            tmp = hc
            for i in xrange(len(pos_poly) - 1):
                tmp = tmp[markers_pair_new[i][1] == 'M']
            tmp[markers_pair_new[-1][1] == 'M'] += 1

            n_cov += 1



    ## Plot minor allele freqs
    #import matplotlib.pyplot as plt
    #plt.figure()
    #plt.scatter(np.arange(len(af2)), af2, s=40, c='k', alpha=0.5)
    #plt.xlabel('Position [bp]')
    #plt.ylabel('Allele freq')
    #plt.xlim(-1, i)
    #plt.ylim(1e-5, 1.05)
    #plt.yscale('log')

    #plt.ion()
    #plt.show()

