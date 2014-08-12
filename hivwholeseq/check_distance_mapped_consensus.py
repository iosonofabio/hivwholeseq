# vim: fdm=indent
'''
author:     Fabio Zanini
date:       04/08/14
content:    Check the distance of mapped reads to their consesus, to spot contamination.
'''
# Modules
import os
import argparse
from operator import itemgetter
import pysam
import numpy as np
from Bio import SeqIO
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.samples import load_sequencing_run, SampleSeq
from hivwholeseq.filenames import get_consensus_filename, get_mapped_filename, \
        get_filter_mapped_summary_filename, get_mapped_suspicious_filename
from hivwholeseq.mapping_utils import get_ind_good_cigars, convert_sam_to_bam,\
        pair_generator, get_range_good_cigars



# Functions
def get_distance_from_reference(ref, reads, pairs=True, threshold=None, VERBOSE=0):
    '''Get the number of mismatches from consensus
    
    Parameters:
       - ref: numpy array of the reference
       - reads: iterator over reads (or pairs) to analyze
       - pairs: analyze in pairs
       - threshold: consider only alleles with quality >= x
    
    Note: deletions and insertions both count 1
    '''
    ds = []
    for reads_i in reads:
        # Get the sum of the distance of the two reads
        d = 0
        if not pairs:
            reads_i = [reads_i]
        for read in reads_i:
            pos_ref = read.pos
            pos_read = 0
            for (bt, bl) in read.cigar:
                if bt == 1:
                    d += 1
                    pos_read += bl
                elif bt == 2:
                    d += 1
                    pos_ref += bl
                elif bt == 0:
                    seqb = np.fromstring(read.seq[pos_read: pos_read + bl], 'S1')
                    diffs = (seqb != ref[pos_ref: pos_ref + bl])
                    if threshold is not None:
                        qualb = np.fromstring(read.qual[pos_read: pos_read + bl], np.int8) - 33
                        diffquals = qualb >= threshold
                        d += (diffs & diffquals).sum()
                    else:
                        d += diffs.sum()
                    pos_ref += bl
                    pos_read += bl
        ds.append(d)
    return np.array(ds, int)


def get_distance_histogram(data_folder, adaID, fragment, maxreads=1000, VERBOSE=0):
    '''Get the distance of reads from their consensus'''
    reffilename = get_consensus_filename(data_folder, adaID, fragment)
    refseq = SeqIO.read(reffilename, 'fasta')
    ref = np.array(refseq)

    bamfilename = get_mapped_filename(data_folder, adaID, fragment, type='bam',
                                      filtered=False)

    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        n_pairs = 0
        read_pairs = []
        for (i, rp) in enumerate(pair_generator(bamfile)):
            if n_pairs >= maxreads:
                break

            r1 = rp[0]
            if not r1.is_proper_pair:
                continue

            read_pairs.append(rp)
            n_pairs += 1

        ds = get_distance_from_reference(ref, read_pairs, threshold=30)

    h = np.bincount(ds)
    return h


def plot_distance_histogram(h, cumulative=True, title='', ax=None, **kwargs):
    '''Plot distance histogram'''
    if ax is None:
        fig, ax = plt.subplots()

    if cumulative:
        x = np.arange(len(h) + 1)
        y = 1.0 - np.insert(1.0 * h.cumsum() / h.sum(), 0, 0)
        y[-1] = min(0.1 * y[-2], 1e-5)
        ax.plot(x, y, lw=2, **kwargs)
        ax.set_ylabel('fraction of read pairs with d > x')
        ax.set_yscale('log')
    else:
        ax.plot(np.arange(len(h)), h, lw=2, **kwargs)
        ax.set_ylabel('# read pairs')

    ax.set_xlabel('Hamming distance from consensus')
    ax.set_xlim(-1, 201)
    ax.grid(True)

    if len(title):
        ax.set_title(title)

    plt.ion()
    plt.show()



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Check distance histogram from consensus',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. 28, 37)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help=('Verbosity level [0-3]'))
    parser.add_argument('--maxreads', type=int, default=1000,
                        help='Maximal number of reads to analyze')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    maxreads = args.maxreads

    # Specify the dataset
    dataset = load_sequencing_run(seq_run)
    data_folder = dataset.folder

    # If the script is called with no adaID, iterate over all
    samples = dataset.samples
    if adaIDs is not None:
        samples = samples.loc[samples.adapter.isin(adaIDs)]
    if VERBOSE >= 3:
        print 'adaIDs', samples.adapter

    hists = []
    for (samplename, sample) in samples.iterrows():
        if VERBOSE == 1:
            print samplename,
        elif VERBOSE >= 2:
            print samplename

        sample = SampleSeq(sample)

        if str(sample.PCR) == 'nan':
            if VERBOSE >= 1:
                print 'PCR type not found, skipping'
            continue

        adaID = sample.adapter

        if not fragments:
            fragments_sample = sample.regions_generic
        else:
            fragments_sample = [fr for fr in fragments if fr in sample.regions_generic]
        if VERBOSE >= 3:
            print 'adaID '+adaID+': fragments '+' '.join(fragments_sample)

        for fragment in fragments_sample:

            bamfilename = get_mapped_filename(data_folder, adaID, fragment, type='bam',
                                              filtered=False)
            if not os.path.isfile(bamfilename):
                if VERBOSE >= 1:
                    print 'missing mapped file, skipping'
                continue

            dist_hist = get_distance_histogram(data_folder, adaID, fragment, VERBOSE=VERBOSE, maxreads=maxreads)
            label = [seq_run, adaID, samplename, fragment, dist_hist.sum()]
            hists.append((dist_hist, label))
    
    if len(hists) == 1:
        label = ', '.join(hists[0][1][:-1])+' ('+str(hists[0][1][-1])+')'
        plot_distance_histogram(hists[0][0], title=label)
    else:
        fig, ax = plt.subplots()
        for i, h in enumerate(hists):
            color = cm.jet(1.0 * i / len(hists))
            label = h[1]
            if len(fragments_sample) == 1:
                del label[3]
            elif len(samples) == 1:
                del label[2]
                del label[1]
            del label[0]
            label = ', '.join(label[:-1])+' ('+str(label[-1])+')'
            plot_distance_histogram(h[0], label=label, ax=ax, color=color)

        ax.legend(loc=1, fontsize=10)
        title = seq_run
        if len(fragments_sample) == 1:
            title = title+', '+fragment
        elif len(samples) == 1:
            title = title+', '+adaID+', '+samplename
        ax.set_title(title)
