# vim: fdm=marker
'''
author:     Fabio Zanini
date:       16/10/13
content:    Find out what the spikes in allele frequency are.
'''
# Modules
import os
import sys
import argparse
import pysam
from collections import defaultdict
from itertools import izip
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt

from mapping.datasets import MiSeq_runs
from mapping.miseq import alpha
from mapping.filenames import get_NL43_entire, get_NL43_fragmented, \
        get_F10_entire, get_F10_fragmented, \
        get_consensus_filename, get_allele_counts_filename, get_coverage_filename, \
        get_mapped_filename
from mapping.mapping_utils import align_muscle, sort_bam, index_bam

from errors_pure_plasmid import references, colors, figures_folder, get_count, \
        get_consensus, get_minor_counts, get_minor_nus, get_coverage, spikes_motifs



# Functions
def fish_reads(miseq_run, adaID, fragment, spike, VERBOSE=0):
    '''Fish out reads from a spike (remember to close the file!)'''
    pos_spike = spike[0]

    # Get BAM file
    data_folder = MiSeq_runs[miseq_run]['folder']
    bamfilename = get_mapped_filename(data_folder, adaID, fragment,
                                      filtered=True, sort=True)
    # Sort BAM if needed
    if not os.path.isfile(bamfilename):
        sort_bam(bamfilename)
        index_bam(bamfilename)

    bamfile = pysam.Samfile(bamfilename, 'rb')
    print bamfile.references
    reference = '_'.join(['adaID', '{:02d}'.format(adaID),
                      fragment, 'consensus'])

    return bamfile, bamfile.fetch(reference, pos_spike, pos_spike+1) 


def troubleshoot_spike(miseq_run, adaID, fragment, spike, reads=None, key=None,
                       VERBOSE=0):
    '''Find the source of high-frequency errors'''
    pos_spike = spike[0]

    if reads is None:
        bamfile, reads = fish_reads(miseq_run, adaID, fragment, spike,
                                    VERBOSE=VERBOSE)

    # Get consensus
    data_folder = MiSeq_runs[miseq_run]['folder']
    consensus =  SeqIO.read(get_consensus_filename(data_folder, adaID, fragment,
                                                   trim_primers=True), 'fasta')
    cons = np.array(consensus)

    # Get allele counts and coverage for testing
    counts = np.load(get_allele_counts_filename(data_folder, adaID, fragment))
    count = get_count(counts)
    cov = get_coverage(count)

    al_WT = consensus[pos_spike]

    ## 1. make sure there are errors
    #counts_pos = np.zeros_like(alpha, int)
    #for read in reads:
    #    # if only fwd/rev are requested, do so
    #    if (key is not None) and (read.is_reverse != (key == 'rev')): continue

    #    # reads are filtered, hence we can proceed to CIGARs
    #    pos = read.pos
    #    seq = read.seq
    #    cl = len(read.cigar)
    #    for ic, (bt, bl) in enumerate(read.cigar):
    #        if bt == 0:
    #            if pos + bl > pos_spike:
    #                al = seq[pos_spike - pos]
    #                counts_pos[(al == alpha).nonzero()[0][0]] += 1
    #                break
    #            if ic != cl - 1:
    #                pos += bl
    #                seq = seq[bl:]
    #        
    #        elif bt == 1:
    #            if ic != cl - 1:
    #                seq = seq[bl:]

    #        elif bt == 2:
    #            if pos + bl > pos_spike:
    #                counts_pos[('-' == alpha).nonzero()[0][0]] += 1
    #                break
    #            pos += bl

    #return al_WT, counts_pos
    ## OK, the error frequency is correct.

    # 2. see where the error is in the read
    pos_read = np.zeros((len(alpha), 250), int)
    pos_read_end = np.zeros_like(pos_read)
    read_names = []
    for read in reads:
        # if only fwd/rev are requested, do so
        if (key is not None) and (read.is_reverse != (key == 'rev')): continue

        # reads are filtered, hence we can proceed to CIGARs
        pos = read.pos
        posr = 0
        seq = read.seq
        cl = len(read.cigar)
        for ic, (bt, bl) in enumerate(read.cigar):
            if bt == 0:
                if pos + bl > pos_spike:
                    al = seq[pos_spike - pos]
                    if al != al_WT:
                        if read.isize < 250:
                            read_names.append(read.qname)
                    pos_read[(al == alpha).nonzero()[0][0], \
                             posr + pos_spike - pos] += 1
                    pos_read_end[(al == alpha).nonzero()[0][0], \
                                 len(read.seq) - 1 - (posr + pos_spike - pos)] += 1
                    break
                if ic != cl - 1:
                    pos += bl
                    posr += bl
                    seq = seq[bl:]
            
            elif bt == 1:
                if ic != cl - 1:
                    seq = seq[bl:]
                    posr += bl

            elif bt == 2:
                if pos + bl > pos_spike:
                    pos_read[('-' == alpha).nonzero()[0][0], posr] += 1
                    pos_read_end[('-' == alpha).nonzero()[0][0], \
                                 len(read.seq) - 1 - posr] += 1
                    break
                pos += bl

    # The end of a read is the start for a reverse read
    if key == 'rev':
        tmp = pos_read.copy()
        pos_read = pos_read_end
        pos_read_end = tmp

    ail_min = np.argsort(pos_read.sum(axis=1))[-2]
    al_min = alpha[ail_min]
    print fragment, pos_spike, key, al_WT, al_min, \
            1.0 * pos_read[(alpha == al_min).nonzero()[0][0]].sum() / pos_read.sum(), \
            str(consensus.seq)[pos_spike-3: pos_spike+4], pos_read_end[ail_min].argmax()

    read_names = np.array(read_names, 'S100')
    read_names.dump('test.npy')

    import ipdb; ipdb.set_trace()
    return al_WT, (pos_read, pos_read_end)
                
    ## 3. Do those reads contain many mutations?
    #n_muts = []
    #pos_muts = defaultdict(int)
    #for read in reads:
    #    # if only fwd/rev are requested, do so
    #    if (key is not None) and (read.is_reverse != (key == 'rev')): continue

    #    # reads are filtered, hence we can proceed to CIGARs
    #    n_mut = 0
    #    pos_mut = set()
    #    pos = read.pos
    #    posr = 0
    #    seq = read.seq
    #    cl = len(read.cigar)
    #    for ic, (bt, bl) in enumerate(read.cigar):
    #        if bt == 0:
    #            if pos + bl > pos_spike:
    #                al = seq[pos_spike - pos]
    #                if al == al_WT:
    #                    n_mut = -1
    #                    break

    #            consloc = cons[pos: pos + bl]
    #            seqb = np.array(list(seq[:bl]))
    #            pos_mut_block = (seqb != consloc).nonzero()[0] + pos
    #            n_mut = len(pos_mut_block)
    #            pos_mut |= set(pos_mut_block)

    #            if ic != cl - 1:
    #                pos += bl
    #                posr += bl
    #                seq = seq[bl:]
    #        
    #        elif bt == 1:
    #            n_mut += 1
    #            pos_mut.add(pos)
    #            if ic != cl - 1:
    #                seq = seq[bl:]
    #                posr += bl

    #        elif bt == 2:
    #            n_mut += bl
    #            pos_mut |= 
    #            pos += bl

    #    if n_mut >= 0:
    #        n_muts.append(n_mut)
    #        for pm in pos_mut:
    #            pos_muts[pm] += 1

    #print fragment, pos_spike, key, al_WT
    #import ipdb; ipdb.set_trace()
    #return al_WT, n_muts



# Script
if __name__ == '__main__':

    # Parse input args: this is used to call itself recursively
    parser = argparse.ArgumentParser(description='Errors in HIV pure strains')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

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
        for fragment in fragments:
            spikes = spikes_motifs(miseq_run, adaID, fragment,
                                   VERBOSE=VERBOSE, plot=False, savefig=False)

            for i in xrange(10):
                for key, spikest in spikes.iteritems():
                    if len(spikest) < i + 1:
                        print 'No spikes ;-)'; continue
                    spike = spikest[i]
                    bamfile, reads = fish_reads(miseq_run, adaID, fragment, spike,
                                                VERBOSE=VERBOSE)

                    al_WT, counts_pos = troubleshoot_spike(miseq_run, adaID,
                                                           fragment, spike,
                                                           reads=reads,
                                                           key=key,
                                                           VERBOSE=VERBOSE)

                    bamfile.close()
