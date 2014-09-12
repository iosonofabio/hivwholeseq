# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/09/14
content:    Build consensus of a patient sample, from the filtered reads.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import ambiguous_dna

from hivwholeseq.patients.patients import load_samples_sequenced, SamplePat
from hivwholeseq.patients.filenames import get_initial_reference_filename



# Function
def pileup_trim_reads_coverfull(bamfile, edges, VERBOSE=0):
    '''Collect reads that fully cover a region, and trim them to the same
    
    Note: this function does not look at the full read pair!
    '''
    pos_ref = edges[0]
    block_len = edges[1] - edges[0]

    seqs = []
    for read in bamfile:
        start_read = end_read = read.pos
        if start_read > pos_ref:
            continue

        end_read += sum(bl for bt, bl in read.cigar if bt in (0, 2))
        if end_read < edges[1]:
            continue

        # Find start of the block in the read
        pos_reft = read.pos
        start_found = False
        pos_read_start = 0
        pos_read_end = 0
        for (bt, bl) in read.cigar:
            if bt == 1:
                if not start_found:
                    pos_read_start += bl
                pos_read_end += bl
            elif bt == 2:
                if (not start_found) and (pos_reft + bl > pos_ref):
                    start_found = True
                if pos_reft + bl > pos_ref + block_len:
                    break
                pos_reft += bl
            else:
                if (not start_found) and (pos_reft + bl > pos_ref):
                    pos_read_start += pos_ref - pos_reft
                    start_found = True
                if pos_reft + bl > pos_ref + block_len:
                    pos_read_end += pos_ref + block_len - pos_reft
                    break

                if not start_found:
                    pos_read_start += bl
                pos_read_end += bl
                pos_reft += bl
        
        seq = read.seq[pos_read_start: pos_read_end]
        seqs.append(seq)

    bamfile.reset()

    return seqs


def pileup_trim_reads_coverstart(bamfile, start, VERBOSE=0):
    '''Collect reads that cover the start of a region, and trim them to the same
    
    Note: this function does not look at the full read pair!
    '''
    pos_ref = start

    seqs = []
    for read in bamfile:
        start_read = end_read = read.pos
        if (start_read > pos_ref) or (start_read < pos_ref - 300):
            continue

        end_read += sum(bl for bt, bl in read.cigar if bt in (0, 2))
        if end_read <= start:
            continue

        # Find start of the block in the read
        pos_reft = read.pos
        pos_read_start = 0
        for (bt, bl) in read.cigar:
            if bt == 1:
                pos_read_start += bl

            elif bt == 2:
                if pos_reft + bl > pos_ref:
                    seq = read.seq[pos_read_start:]
                    break

                pos_reft += bl
            else:
                if pos_reft + bl > pos_ref:
                    pos_read_start += pos_ref - pos_reft
                    seq = read.seq[pos_read_start:]
                    break

                pos_read_start += bl
                pos_reft += bl
        
        seqs.append(seq)

    bamfile.reset()

    return seqs


def join_block_to_consensus(consensus, cons_block, seed_length=40, VERBOSE=0):
    '''Join a new block to an extant consensus'''
    import numpy as np

    start = consensus.rfind(cons_block[:seed_length])
    if start == -1:
        # Rescue for single polymorphisms
        start_tail = max(0, len(consensus) - 300)
        cons_tail = np.fromstring(consensus[start_tail:], 'S1')
        seed = np.fromstring(cons_block[:seed_length], 'S1')
        sl = len(seed)
        n_matches = [(seed == cons_tail[i: i + sl]).sum() for i in xrange(len(cons_tail) - sl)]
        pos = np.argmax(n_matches)
        if n_matches[pos] < 0.80 * sl:

            # Try a local pairwise alignment
            seed = ''.join(seed)
            cons_tail = ''.join(cons_tail)
            from seqanpy import align_local
            scores = []
            seeds1 = []
            seeds2 = []
            for i in xrange(len(cons_tail) - sl):
                score, seed1, seed2 = align_local(cons_tail[i: i + sl], seed)
                scores.append(score)
                seeds1.append(seed1)
                seeds2.append(seed2)

            ind = np.argmax(scores)
            scoremax = scores[ind]
            if scoremax < 20:
                raise ValueError('Hole in consensus')
            
            pos = cons_tail.rfind(seeds1[ind].replace('-', ''))
            cons_block = cons_block[cons_block.find(seeds2[ind].replace('-', '')):]

        start = start_tail + pos

    consensus = consensus[:start] + cons_block

    return consensus


def build_consensus(bamfilename, len_reference, VERBOSE=0,
                    block_len=100,
                    reads_per_alignment=31):
    '''Build a consensus from mapped filtered reads'''
    if VERBOSE:
        print 'Build consensus'
    
    from operator import itemgetter
    import numpy as np
    import pysam
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet.IUPAC import ambiguous_dna

    from hivwholeseq.miseq import alpha
    from hivwholeseq.mapping_utils import pair_generator
    from hivwholeseq.sequence_utils import build_local_consensus
    
    with pysam.Samfile(bamfilename, 'rb') as bamfile:

        # Get first block covered, even if partially, and record where each read started
        block_len = block_len
        seqs = []
        n_block = 0
        while not seqs:
            start_block = n_block * (block_len // 2)
            for read in bamfile:
                if read.pos <= start_block:
                    seqs.append((read.pos, ('N' * read.pos) + read.seq[:block_len - read.pos]))
            bamfile.reset()
            n_block += 1
        
        # If there are too many reads, take the reads that start earliest
        if len(seqs) > reads_per_alignment:
            np.random.shuffle(seqs)
            seqs.sort(key=itemgetter(0))
            seqs = seqs[:reads_per_alignment]

        seqrecs = [SeqRecord(Seq(s, ambiguous_dna), id=str(i), name=str(i))
                   for i, (pos, s) in enumerate(seqs)]
        consensus = build_local_consensus(seqrecs, VERBOSE=VERBOSE, full_cover=False)

        # Block, by block, make local alignment and join to previous consensus
        # There are two ways of finishing the loop:
        # 1. if we cover all the way to the end of the reference, good
        # 2. if we find no reads fully covering a block BEFORE that, add a final block
        while start_block < len_reference:
            edges = (start_block, min(len_reference, start_block + block_len))

            if VERBOSE >= 2:
                print 'block n.', n_block, 'region:', edges

            seqs = pileup_trim_reads_coverfull(bamfile, edges, VERBOSE=VERBOSE)

            # If we do not find reads that fully cover, consider it the end of
            # the consensus, only the final block is missing
            if not seqs:
                break
            elif len(seqs) > reads_per_alignment:
                np.random.shuffle(seqs)
                seqs = seqs[:reads_per_alignment]

            # Local alignment
            # --------------
            # -----   ------
            # --------   ---
            #---------------
            seqrecs = [SeqRecord(Seq(s, ambiguous_dna), id=str(i), name=str(i))
                       for i, s in enumerate(seqs)]
            cons_block = build_local_consensus(seqrecs, VERBOSE=VERBOSE, full_cover=True)

            # Join to the rest of the consensus, like this:
            # ---------------------------
            #                          ---------
            consensus = join_block_to_consensus(consensus, cons_block,
                                                seed_length=block_len // 4,
                                                VERBOSE=VERBOSE)

            start_block += 2 * block_len // 3
            n_block += 1

        # If we cover the whole reference, good
        else:
            return consensus

        if VERBOSE >= 2:
            print 'final block'

        # If we broke out of the while, a final block is needed
        seqs = pileup_trim_reads_coverstart(bamfile, start_block, VERBOSE=VERBOSE)

        # Sort reads by length
        if len(seqs) > reads_per_alignment:
            np.random.shuffle(seqs)
            seqs.sort(key=len, reverse=True)
            seqs = seqs[:reads_per_alignment]

            # Complete with N, approximately
            sl = len(seqs[0])
            seqs = [s+('N' * (sl - len(s))) for s in seqs]


        seqrecs = [SeqRecord(Seq(s, ambiguous_dna), id=str(i), name=str(i))
                   for i, s in enumerate(seqs)]
        cons_block = build_local_consensus(seqrecs, VERBOSE=VERBOSE, full_cover=False)
        consensus = join_block_to_consensus(consensus, cons_block, VERBOSE=VERBOSE)

    return consensus



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Build consensus sequence for a patient sample',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=True)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to map (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save the allele counts to file')
    parser.add_argument('--block-length', type=int, default=200, dest='block_len',
                        help='Length of each local consensus block')
    parser.add_argument('--reads-per-alignment', type=int, default=31,
                        dest='reads_per_alignment',
                        help='Number of (random) reads used for the local consensi')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    fragments = args.fragments
    VERBOSE = args.verbose
    save_to_file = args.save
    n_reads_per_ali = args.reads_per_alignment
    block_len = args.block_len

    samples = load_samples_sequenced()
    if pnames is not None:
        samples = samples.loc[samples.patient.isin(pnames)]
    elif samplenames is not None:
        samples = samples.loc[samples.index.isin(samplenames)]

    if VERBOSE >= 2:
        print 'samples', samples.index.tolist()

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    counts_all = []
    for fragment in fragments:
        counts = []
        for samplename, sample in samples.iterrows():
            if VERBOSE >= 1:
                print samplename, fragment,
                if VERBOSE >= 2:
                    print ''

            sample = SamplePat(sample)
            pname = sample.patient
            refseq = SeqIO.read(get_initial_reference_filename(pname, fragment), 'fasta')
            refm = np.array(refseq)
            len_reference = len(refseq)

            for PCR in (1, 2):
                bamfilename = sample.get_mapped_filtered_filename(fragment, PCR=PCR)
                if not os.path.isfile(bamfilename):
                    continue
                
                if VERBOSE >= 1:
                    print 'PCR', PCR,
                    if VERBOSE >= 2:
                        print ''

                cons = build_consensus(bamfilename, len_reference, VERBOSE=VERBOSE,
                                       block_len=block_len,
                                       reads_per_alignment=n_reads_per_ali)
                consm = np.fromstring(cons, 'S1')

                if VERBOSE >= 2:
                    print 'Reference length:', len_reference, 'consensus:', len(cons),
                    if len_reference != len(cons):
                        from seqanpy import align_local
                        (score, ali1, ali2) = align_local(''.join(refseq), cons)
                        alim1 = np.fromstring(ali1, 'S1')
                        alim2 = np.fromstring(ali2, 'S1')
                        n_diff = (alim1 != alim2).sum()
                        print 'overlap', len(alim1), 'n diffs.', n_diff
                    else:
                        n_diff = (refm != consm).sum()
                        print 'n diffs.', n_diff

                if save_to_file:
                    if VERBOSE >= 2:
                        print 'Save to file'

                    fn_out = sample.get_consensus_filename(fragment, PCR=PCR)
                    consrec = SeqRecord(Seq(cons, ambiguous_dna),
                                        id=samplename+'_consensus',
                                        name=samplename+'_consensus',
                                        description=samplename+', consensus',
                                       )
                    SeqIO.write(consrec, fn_out, 'fasta')

            if VERBOSE == 1:
                print ''
