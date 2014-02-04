#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       29/01/14
content:    Get (and plot) the allele counts and frequencies.
'''
# Modules
import os
import argparse
from itertools import izip
from Bio import SeqIO
import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from hivwholeseq_pacbio.samples import samples
from hivwholeseq_pacbio.filenames import get_premapped_file, \
        get_reference_premap_filename


# Globals
data_folder_dict = {'Upp23': '/ebio/ag-neher/share/data/PacBio_HIV_Karolinska/run23/'}
alpha = np.array(['A', 'C', 'G', 'T', 'N', '-'])
alphal = list(alpha)



# Functions
def get_allele_counts(bamfilename, reflen, qual_min=20, maxreads=-1, VERBOSE=0):
    '''Get the allele counts from the reads'''

    counts = np.zeros((len(alpha), reflen), int)
    #TODO: insertions

    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        for i, read in enumerate(bamfile):
            if VERBOSE >= 2:
                if not ((i+1) % 100):
                    print (i+1)

            if i == maxreads:
                break

            if read.is_unmapped:
                if VERBOSE >= 3:
                    print 'Unmapped'
                continue

            # Skip obvious mismappings
            if sum(bl for (bt, bl) in read.cigar if bt in (1, 2)) > 100:
                if VERBOSE >= 3:
                    print 'Many insertions/deletions'
                continue

            seq = np.fromstring(read.seq, 'S1')
            qual = np.fromstring(read.qual, np.int8) - 33
            pos_ref = read.pos
            pos_read = 0

            for (bt, bl) in read.cigar:

                # Matches
                if bt == 0:
                    seqb = seq[pos_read: pos_read + bl]
                    qualb = qual[pos_read: pos_read + bl]
        
                    # Increment counts
                    for j, a in enumerate(alpha):
                        posa = ((seqb == a) & (qualb >= qual_min)).nonzero()[0]
                        if len(posa):
                            counts[j, pos_ref + posa] += 1
        
                    pos_ref += bl
                    pos_read += bl

                # Deletions (have no quality! - FIXME: looks at neighbours?)
                elif bt == 2:
                    counts[-1, pos_ref + np.arange(bl)] += 1
                    pos_ref += bl

                # TODO: insertions
                elif bt == 1:
                    pos_read += bl

    return counts


def get_consensus(counts):
    consm = np.zeros(counts.shape[-1], 'S1')
    for i, count in enumerate(counts.T):
        consm[i] = alpha[np.argmax(count)]

    return consm.tostring()


def get_minor_nu(counts, cov_min=100):
    '''Get minor allele frequencies from counts'''
    num = np.ma.masked_all(counts.shape[-1])
    for i, count in enumerate(counts.T):
        # Exclude deletions and see
        count = count[:-1]
        if np.ma.is_masked(count) or (count.sum() < 100):
            continue

        nu = 1.0 * np.sort(count)[-2] / count.sum()
        num[i] = nu
    return num


def plot_minor_frequencies(counts, title):
    '''Plot the coverage'''

    fig, ax = plt.subplots(1, 1, figsize=(14, 6))

    num = get_minor_nu(counts)

    ax.plot(num, color='k', lw=1.5)
    ax.scatter(np.arange(counts.shape[-1]),
               num,
               s=40,
               edgecolor='none', color='k', lw=1.5)
    ax.set_xlabel('Position [bp]')
    ax.set_ylabel('Minor allele freq')

    ax.set_xlim(-50, counts.shape[-1] + 50) 
    ax.set_ylim(1e-5, 1)
    ax.set_yscale('log')
    ax.set_title(title)

    plt.ion()
    plt.show()


def annotate_plot_PCR_primers(ax, refseq, VERBOSE=0):
    '''Annotate the plot with the PCR primer positions'''
    from matplotlib.patches import Rectangle
    from matplotlib import cm
    from hivwholeseq.primer_info import primers_PCR
    from hivwholeseq.annotate_genomewide_consensus import annotate_sequence as aseq

    # Get axes limits
    (ymin, ymax) = ax.get_ylim()
    if ax.get_yscale() == 'linear':
        yspan = ymax - ymin
        yannos = [ymin - (0.05 + 0.03 * i) * yspan for i in xrange(len(primers_PCR))]
        ymin_new = ymin - 0.4 * yspan
    else:
        yspan = ymax / ymin
        yannos = [np.exp(np.log(ymin) - (0.05 + 0.03 * i) * np.log(yspan))
                  for i in xrange(len(primers_PCR))]
        ymin_new = ymin / np.exp(0.4 * np.log(yspan))

    # Plot annotations
    aseq(refseq, features=['PCR primers'])
    for (name, (pr_fwd, pr_rev)) in primers_PCR.iteritems():
        # Only F5a used
        if 'F5a' in name:
            continue

        feat = [fea for fea in refseq.features if fea.id == name][0]
        fwd_start = feat.location.nofuzzy_start

        # backwards
        rev_start = feat.location.nofuzzy_end - len(pr_rev)

        # Annotate fwd primer
        ax.add_patch(Rectangle((fwd_start, ymin_new), len(pr_fwd), ymax - ymin_new,
                               fc='b',
                               alpha=0.5,
                               ec='none'))
        ax.plot([fwd_start, fwd_start + len(pr_fwd)], 2 * [yannos[i]], lw=2.5,
                 c='k', alpha=0.5)
        ax.text(fwd_start + 0.1 * len(pr_fwd), yannos[i], feat.id)

        # Annotate rev primer
        ax.add_patch(Rectangle((rev_start, ymin_new), len(pr_rev), ymax - ymin_new,
                               fc='r',
                               alpha=0.5,
                               ec='none'))
        ax.plot([rev_start, rev_start + len(pr_rev)], 2 * [yannos[i]], lw=2.5,
                 c='k', alpha=0.5)
        ax.text(rev_start + 0.1 * len(pr_rev), yannos[i], feat.id)

    plt.draw()


def remove_primers(counts, refseq, VERBOSE=0):
    '''Mask the allele counts at the primer positions'''
    from hivwholeseq.primer_info import primers_PCR
    counts_noprim = np.ma.array(counts, int)
    for feat in refseq.features:
        if 'F5a' in feat.id:
            continue

        (pr_fwd, pr_rev) = primers_PCR[feat.id]

        fwd_start = feat.location.nofuzzy_start
        fwd_end = fwd_start + len(pr_fwd)

        rev_start = feat.location.nofuzzy_end - len(pr_rev)
        rev_end = rev_start + len(pr_rev)

        counts_noprim[:, fwd_start: fwd_end] = np.ma.masked
        counts_noprim[:, rev_start: rev_end] = np.ma.masked

    return counts_noprim



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get PacBio allele counts and \
                                     frequencies')
    parser.add_argument('--run', default='Upp23',
                        help='PacBio run to analyze (e.g. Upp23)')
    parser.add_argument('--sample', required=True,
                        help='Sample to analyze (e.g. S1)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to map')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')

    args = parser.parse_args()
    seq_run = args.run
    samplename = args.sample
    VERBOSE = args.verbose
    maxreads = args.maxreads
    submit = args.submit

    # Specify the dataset
    data_folder = data_folder_dict[seq_run]
    sample = samples.set_index('name').loc[samplename]

    # Get NL4-3 reference
    refseq = SeqIO.read(get_reference_premap_filename(data_folder, samplename), 'fasta')

    # Get counts
    allele_counts_filename = data_folder+samplename+'/allele_counts.npy'
    try:
        counts = np.load(allele_counts_filename)
    except IOError:
        bamfilename = get_premapped_file(data_folder, samplename)
        counts = get_allele_counts(bamfilename, len(refseq), maxreads=maxreads,
                                   VERBOSE=VERBOSE,
                                   qual_min=33)
        counts.dump(allele_counts_filename)

    # Plot
    plot_minor_frequencies(counts + 1.5e-5, 'PacBio minor allele freq: '+samplename)

    annotate_plot_PCR_primers(plt.gca(), refseq, VERBOSE=VERBOSE)

    # Replot without those peaks (masked)
    counts_noprim = remove_primers(counts, refseq, VERBOSE=VERBOSE)
    plot_minor_frequencies(counts_noprim + 1.5e-5, 'PacBio minor allele freq: '+samplename+\
                           ', PCR primers sites removed')


