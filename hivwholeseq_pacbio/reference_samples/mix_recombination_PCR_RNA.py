# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/09/14
content:    Characterize recombination levels in conventional vs emulsion RT-PCR.
'''
# Modules
import sys
import os
import argparse
import numpy as np
from collections import Counter
from operator import itemgetter, attrgetter
from Bio import SeqIO, AlignIO
import pysam
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.mapping_utils import align_muscle, reads_to_seqrecord
from hivwholeseq.generic_utils import getchar
from hivwholeseq_pacbio.sequencing.pacbio_rs_II import alphal
from hivwholeseq_pacbio.sequencing.samples import load_sample
from hivwholeseq_pacbio.reference import load_custom_reference



# Globals
fragment = 'F4'
mix_samples = {'Upp91': {'38540': 0.495, '38304': 0.495, 'LAI-III': 0.01}}


# Functions
def analyze_recombination_read(refseqs, read, VERBOSE=0):
    '''Analyze recombination haplotype of a single PacBio read'''
    seqs = refseqs + reads_to_seqrecord([read])
    ali = align_muscle(*seqs, sort=True)

    # Trim left edge
    if '-' in ali[:-1, 0]:
        start = max(len(seq) - len(''.join(seq).lstrip('-')) for seq in ali[:-1])
    else:
        start = 0

    # Trim right edge
    if '-' in ali[:-1, -1]:
        end = min(len(''.join(seq).rstrip('-')) for seq in ali[:-1])
    else:
        end = ali.get_alignment_length()

    ali = ali[:, start: end]

    # Look for polymorphisms
    alim = np.array(ali)
    read_bitset = np.zeros(alim.shape[1], np.uint8)
    for i, row in enumerate(alim[:-1]):
        read_bitset += ((alim[-1] == row) << i)
    return read_bitset



# Script
if __name__ == '__main__':

    # Parse input args: this is used to call itself recursively
    parser = argparse.ArgumentParser(description='Errors in HIV pure strains')
    parser.add_argument('--sample', required=True,
                        help='Sample to study (e.g. convPCR)')
    parser.add_argument('--maxreads', type=int, default=100,
                        help='Number of reads analyzed')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    samplename = args.sample
    maxreads = args.maxreads
    VERBOSE = args.verbose

    sample = load_sample(samplename)
    seq_run = sample.run
    data_folder = sample.seqrun_folder

    mix_samplenames = mix_samples[seq_run] 
    refseqs = [load_custom_reference(refn+'_'+fragment) for refn in mix_samplenames]

    # Use only two refs for now
    refseqs = filter(lambda x: 'LAI' not in x.id, refseqs)

    print map(attrgetter('name'), refseqs)
    frac_all = Counter()
    hist_ratio = []

    bamfilename = sample.get_premapped_filename()
    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        for ir, read in enumerate(bamfile):
            if ir == maxreads:
                break

            if VERBOSE == 2:
                if not ((ir + 1) % 10):
                    print ir + 1

            read_bitset = analyze_recombination_read(refseqs, read, VERBOSE=VERBOSE) 

            # Check whether it's mapped, strictly
            if (read_bitset != 0).sum() < 1000:
                continue

            # Exclude if all or none have it
            ind_poly = (read_bitset != 0) & (read_bitset + 1 < (1 << len(refseqs)))
            read_bitset_poly = read_bitset[ind_poly]
            frac_read = Counter(read_bitset_poly)

            frac_all += frac_read
            if (len(read_bitset_poly) >= 50):
                # FIXME: this only works for 2 refs only
                ratio = 1.0 * frac_read[1] / sum(frac_read.itervalues())
                hist_ratio.append(ratio)

            if VERBOSE >= 3:
                print ir + 1

                if VERBOSE >= 4:
                    print ''.join(map(str, read_bitset))
                else:
                    print ''.join(map(str, read_bitset_poly))

                ch = getchar()
                if ch.lower() == 'q':
                    break

                print ''


    # Print ratio of allele from refseq[0] / refseq[1]
    hist_ratio_cum = np.sort(hist_ratio)
    fig, ax = plt.subplots()
    ax.plot(hist_ratio_cum, 1.0 - np.linspace(0, 1, len(hist_ratio_cum)), lw=2)
    ax.set_xlabel('Allele ratio in read: '+refseqs[0].name.split('_')[0]+\
                  ' / '+refseqs[1].name.split('_')[0])
    ax.set_ylabel('fraction of reads with ratio > x')
    ax.set_title(seq_run+', '+samplename)
    plt.tight_layout()

    plt.ion()
    plt.show()
