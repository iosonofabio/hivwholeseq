# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/10/14
content:    Build a local phylogenetic tree from the reads.
'''
# Modules
import argparse
import numpy as np
from Bio import Phylo

from hivwholeseq.mapping_utils import align_muscle
from hivwholeseq.patients.patients import load_patient
from hivwholeseq.argparse_utils import RoiAction



# Functions
def get_local_blocks(reads, start, end, VERBOSE=0):
    '''Extract read pairs covering a genomic region, possibly trimming them to it'''
    from hivwholeseq.mapping_utils import get_read_start_end

    pairs_out = []
    for irp, read_pair in enumerate(reads):
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
        coords = map(get_read_start_end, reads)
        # Either a single read covers the whole thing...
        cover = np.array([False, False], bool)
        for ir, (start_read, end_read) in enumerate(coords):
            if (start_read <= start) and (end_read >= end):
                cover[ir] = True

        # or require 20 nucleotides overlap, to be able to align the two reads





        if start_fwd > start:
            continue
        if end_rev < end:
            continue
        if (end_fwd < end) and (start_rev > start) and (start_rev > end_fwd):
            continue
        # FIXME: require a certain overlap if both cover
        # Special case: if they exactly lie next ot each other but with insertions,
        # it's a big mess, discard for now
        if (end_fwd == start_rev) and \
           ((reads[0].cigar[-1][0] == 1) or (reads[1].cigar[0][0] == 1)):
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
            for (bt, bl) in read.cigar:
                if bt == 1:
                    if pos_ref >= start:
                        seq.append(read.seq[pos_read: pos_read + bl])
                    pos_read += bl

                elif bt == 2:
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
            seqs.append(seq)

        # Merge sequences if both fwd and rev cover part of it
        if seqs[0] is None:
            seq = seqs[1]
        elif seqs[1] is None:
            seq = seqs[0]
        else:
            # The fwd read starts before the rev, because of our insert sizes
            # Align them and merge
            if end_fwd == start_rev:
                seq = seqs[0] + seqs[1]
            elif (end_fwd - start_rev) < 20:


            else:
                seq1 = np.fromstring(seqs[0], 'S1')
                seed = seqs[1][:20]
                sl = len(seed)
                n_matches = [(seed == seq1[i: i + sl]) for i in xrange(len(seq1) - sl)]
                pos = len(n_matches) - 1 - np.argmax(n_matches[::-1])
                if n_matches[pos] < 0.75 * sl:
                    continue




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




# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get local tree from reads',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--roi', required=True, action=RoiAction,
                        help='Region of interest (e.g. F1 300 350)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--PCR1', type=int, default=1,
                        help='Take only PCR1 samples [0=off, 1=where both available, 2=always]')
    parser.add_argument('--maxreads', type=int, default=10,
                        help='Number of reads analyzed per sample')
    parser.add_argument('--plot', action='store_true',
                        help='Plot tree')

    args = parser.parse_args()
    pname = args.patient
    roi = args.roi
    VERBOSE = args.verbose
    use_PCR1 = args.PCR1
    maxreads = args.maxreads
    use_plot = args.plot
