# vim: fdm=marker
'''
author:     Fabio Zanini
date:       04/09/14
content:    Get local haplotypes from the reads fully covering a small region
            of interest, cluster them somewhat and pack them into trajectories.

            NOTE: this script ignores insertions, which require a multiple
            sequence alignment, for the sake of simplicity. Large regions are
            going to fail, because read pairs cover them with a hole in the middle.
'''
# Modules
import os
import argparse
from operator import itemgetter
from itertools import izip
import numpy as np
from Bio import SeqIO

from hivwholeseq.patients.patients import load_patient



# Functions
def get_local_block(bamfilename, start, end, VERBOSE=0, maxreads=-1, refroi=None):
    '''Extract reads fully covering the region, discarding insertions'''
    import sys
    import pysam
    from hivwholeseq.mapping_utils import pair_generator
    from hivwholeseq.mapping_utils import extract_mapped_reads_subsample_open

    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        block = []

        if maxreads == -1:
            reads_iter = pair_generator(bamfile)
        else:
            reads_iter =  extract_mapped_reads_subsample_open(bamfile, maxreads,
                                                              VERBOSE=VERBOSE,
                                                              pairs=True)

        for irp, reads in enumerate(reads_iter):
            if VERBOSE >= 2:
                if not ((irp + 1) % 10000):
                    if irp + 1 != 10000:
                        sys.stdout.write("\x1b[1A\n")
                    sys.stdout.write(str(irp + 1))
                    sys.stdout.flush()

            # Sort fwd read first
            is_fwd = reads[0].is_reverse
            reads = [reads[is_fwd], reads[not is_fwd]]

            # Check for coverage of the region
            start_fwd = reads[0].pos
            end_fwd = start_fwd + sum(bl for (bt, bl) in reads[0].cigar if bt in (0, 2))
            start_rev = reads[1].pos
            end_rev = start_rev + sum(bl for (bt, bl) in reads[1].cigar if bt in (0, 2))
            if start_fwd > start:
                continue
            if end_rev < end:
                continue
            if (end_fwd < end) and (start_rev > start) and (start_rev > end_fwd):
                continue

            if VERBOSE >= 3:
                print ' '.join(map('{:>4d}'.format,
                                   [start_fwd, end_fwd, start_rev, end_rev]))

            # Gather info from both reads, merge by putting ambiguous nucleotides
            seqs = []
            st_ens = [[start_fwd, end_fwd], [start_rev, end_rev]]
            for ir, read in enumerate(reads):
                (start_read, end_read) = st_ens[ir]
                if (end_read <= start) or (start_read >= end):
                    seqs.append(None)
                    continue

                seq = []
                pos_ref = start_read
                pos_read = 0
                start_block = max(start, start_read) - start
                end_block = min(end, end_read) - start
                for (bt, bl) in read.cigar:
                    if bt == 1:
                        pos_read += bl

                    elif bt == 2:
                        if pos_ref + bl > start:
                            st = max(0, start - pos_ref)
                            en = min(bl, end - pos_ref)
                            seq.append('-' * (en - st))
                            if pos_ref + bl >= end:
                                break
                        pos_ref += bl

                    elif bt == 0:
                        if pos_ref + bl > start:
                            st = max(0, start - pos_ref)
                            en = min(bl, end - pos_ref)
                            seq.append(read.seq[pos_read + st: pos_read + en])
                            if pos_ref + bl >= end:
                                break
                        pos_ref += bl
                        pos_read += bl
                seq = ''.join(seq)
                seqs.append((start_block, end_block, seq))

            # Merge sequences if both fwd and rev cover part of it
            if seqs[0] is None:
                seq = seqs[1][2]
            elif seqs[1] is None:
                seq = seqs[0][2]
            else:
                # The fwd read starts before the rev, because of our insert sizes
                end_block_fwd = seqs[0][1]
                start_block_rev = seqs[1][0]
                overlap = [seqs[0][2][start_block_rev:],
                           seqs[1][2][:end_block_fwd - start_block_rev]]

                # The two reads in a pair should have the same length in the overlap
                if len(overlap[0]) != len(overlap[1]):
                    if VERBOSE >= 3:
                        print 'WARNING:', reads[0].qname, 'not same length in overlap!'
                    continue

                ol_fwd = np.fromstring(overlap[0], 'S1')
                ol_rev = np.fromstring(overlap[1], 'S1')

                ol_fwd[ol_fwd != ol_rev] = 'N'
                seq = seqs[0][2][:start_block_rev] + \
                      ol_fwd.tostring() + \
                      seqs[1][2][end_block_fwd - start_block_rev:]

            block.append(seq)

        if VERBOSE >= 2:
            print ''

    return block


def cluster_block(block, VERBOSE=0):
    '''Cluster sequences in a block'''
    from collections import Counter
    from itertools import product

    # Identical sequences are clustered first
    clusters = Counter(block)
    
    # Cluster ambiguous nucleotides
    # TODO: this is buggy at present, and reads with Ns are a small minority
    # (typically < 2%), without any bias towards sequence, so do we want to do this?
    clusters_new = Counter()
    n_reads_clu = {'no N': 0, 'with N': 0}
    for (seq, count) in clusters.iteritems():
        if 'N' not in seq:
            clusters_new[seq] = count
            n_reads_clu['no N'] += count
        else:
            n_reads_clu['with N'] += count

        # Skip if too many Ns (exponentially many neighbouring seqs)
        #nN = seq.count('N')
        #if nN > 3:
        #    continue

        #alpha = ['A', 'C', 'G', 'T']
        #best_neighbour = None
        #count_neighbour = 0
        #for comb in product(*([alpha] * nN)):
        #    seqtmp = seq
        #    for nuc_new in comb:
        #        seqtmp = seqtmp.replace('N', nuc_new, 1)
        #    count_neighbour_tmp = clusters[seqtmp]
        #    if count_neighbour_tmp > count_neighbour:
        #        best_neighbour = seqtmp
        #        count_neighbour = count_neighbour_tmp

        #if best_neighbour is not None:
        #    clusters_new[best_neighbour] += count
    
    if VERBOSE >= 3:
        print '# clusters: no N', n_reads_clu['no N'], 'with N', n_reads_clu['with N']

    clusters = clusters_new
    
    # Remove zeros (there should be none...)
    clusters_new += Counter()

    return clusters


def most_common_fractions(ali):
    '''Find most common haplotype and its frequency'''
    mcf = []
    for a in ali:
        a = a[2]
        (seq, count) = a.most_common(1)[0]
        frac = 1.0 * count / sum(a.itervalues())
        mcf.append([seq, frac])

    return mcf 



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get coverage trajectories',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--roi', required=True, nargs='+',
                        help='Region of interest (e.g. F1 300 350)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--PCR1', type=int, default=1,
                        help='Take only PCR1 samples [0=off, 1=where both available, 2=always]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of reads analyzed per sample')
    parser.add_argument('--plot', action='store_true',
                        help='Plot local haplotype trajectories')

    args = parser.parse_args()
    pname = args.patient
    roi = args.roi
    VERBOSE = args.verbose
    use_PCR1 = args.PCR1
    maxreads = args.maxreads
    use_plot = args.plot

    if len(roi) % 3:
        raise ValueError('roi syntax: --roi FRAGMENT START END')
    fragment = roi[0]
    start = int(roi[1])
    end = int(roi[2])

    patient = load_patient(pname)
    patient.discard_nonsequenced_samples()
    refseq = patient.get_reference('F4')
    refroi = ''.join(refseq[start: end])

    if VERBOSE >= 1:
        print patient.name, fragment, start, end

    ali = []  
    for t, sample in izip(patient.times, patient.itersamples()):
        if VERBOSE >= 1:
            print t, sample.name

        bamfilename1 = sample.get_mapped_filtered_filename(fragment, PCR=1)
        bamfilename2 = sample.get_mapped_filtered_filename(fragment, PCR=2)

        if use_PCR1 == 0:
            bamfilenames = []
            if os.path.isfile(bamfilename1):
                bamfilenames.append((1, bamfilename1))
            if os.path.isfile(bamfilename2):
                bamfilenames.append((2, bamfilename1))
        elif use_PCR1 == 1:
            if os.path.isfile(bamfilename1):
                bamfilenames = [(1, bamfilename1)]
            elif os.path.isfile(bamfilename2):
                bamfilenames = [(2, bamfilename2)]
        else:
            if os.path.isfile(bamfilename1):
                bamfilenames = [(1, bamfilename1)]
            else:
                bamfilenames = []
            
        if not bamfilenames:
            continue

        for (PCR, bamfilename) in bamfilenames:
            block = get_local_block(bamfilename, start, end, VERBOSE=VERBOSE,
                                    maxreads=maxreads, refroi=refroi)
            clusters = cluster_block(block)
            ali.append((t, PCR, clusters))

    mcfs = most_common_fractions(ali)
    print 'Most abundant haplotypes'
    print '{:^6}'.format('time'), ('{:^'+str(end - start)+'}').format('haplotype'), 'fraction'
    print '-' * (6 + 1 + end - start + 1 + 8)
    for (t, mcf) in izip(patient.times, mcfs):
        print '{:>6.1f}'.format(t),
        print mcf[0],
        print '{:>8.0%}'.format(mcf[1])


    if use_plot:
        haplos = set()
        for (t, PCR, clusters) in ali:
            haplos |= set(clusters)
        haplos = list(haplos)
        haplos = np.array(haplos, 'S'+str(len(haplos[0])))

        ts = []
        hct = np.zeros((len(ali), len(haplos)), int)
        for it, (t, PCR, clusters) in enumerate(ali):
            ts.append(t)
            for ih, haplo in enumerate(haplos):
                hct[it, ih] = clusters[haplo]

        hft = np.zeros_like(hct, float)
        for it, hc in enumerate(hct):
            hft[it] = 1.0 * hc / hc.sum()

        # Sort by size
        ind = hft.sum(axis=0).argsort()[::-1]
        haplos = haplos[ind]
        hct = hct[:, ind]
        hft = hft[:, ind]

        # Take only major ones and compress the rest
        n_haplos = 50
        haplos_comp = haplos[:n_haplos].tolist() + ['other']
        hft_comp = np.zeros((len(ali), n_haplos + 1), float)
        for ihf in xrange(n_haplos):
            hft_comp[:, ihf] = hft[:, ihf]
        
        hft_comp[:, -1] = hft[:, n_haplos:].sum(axis=1)

        from matplotlib import cm
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        colors = [cm.jet(1.0 * ihf / (n_haplos + 1)) for ihf in xrange(n_haplos + 1)]
        ax.stackplot(ts, *(hft_comp.T), baseline='zero',
                     label=haplos_comp,
                     colors=colors)

        ax.set_xlabel('Time from transmission [days]')
        ax.set_ylabel('Haplotype frequency')
        ax.set_title(pname+', roi = '+' '.join(roi))
        ax.legend(loc=1)
        ax.set_ylim(-0.01, 1.01)
        ax.set_xlim(ts[0] - 20, ts[-1] + 20)

        plt.ion()
        plt.show()





