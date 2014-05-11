#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/03/14
content:    Build consensus after premap, trim and divide. Thanks to the premap
            we have an easy time here: we collect reads a bit from all over the
            fragment and make local consensi (100 bp), which we then chain.

            This script is not working in case of unexpected insertions/deletions,
            e.g. the 1kb+ plasmid insertion in SF162. De novo assembly is needed
            for that, which we could script but for now looks quite useless: the
            divided reads are remapped to initial patient references anyway.
'''
# Modules
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio import AlignIO

from hivwholeseq.mapping_utils import align_muscle
from hivwholeseq.samples import samples
from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.filenames import get_divided_filename, \
        get_reference_premap_filename, \
        get_consensus_filename, \
        get_allele_counts_filename, \
        get_build_consensus_summary_filename
from hivwholeseq.fork_cluster import fork_build_consensus as fork_self



# Functions
def get_reference_all_filename(data_folder, adaID, fragment, ext=True):
    '''Get the file with the cumulated consensi'''
    from hivwholeseq.filenames import foldername_adapter
    fn = '_'.join(['consensus', 'reference', 'ali', fragment])
    fn = data_folder+foldername_adapter(adaID)+fn
    if ext:
        fn = fn+'.fasta'
    return fn


def build_local_consensus(ali, VERBOSE=0, store_allele_counts=False):
    '''Build a local consensus from an MSA'''
    import numpy as np
    from hivwholeseq.miseq import alpha

    allele_counts = np.array([(ali == a).sum(axis=0) for a in alpha], int, ndmin=2)
    cons_local = []
    for counts in allele_counts.T:
        # Pick a random nucleotide in case of a tie
        maxinds = (counts == counts.max()).nonzero()[0]
        maxind = maxinds[np.random.randint(len(maxinds))]
        cons_local.append(alpha[maxind])
    cons_local = np.array(cons_local, 'S1')

    ind_nongap = cons_local != '-'
    cons_local = ''.join(cons_local[ind_nongap])
    
    if store_allele_counts:
        allele_counts = allele_counts[:, ind_nongap]
        return (cons_local, allele_counts)

    return cons_local


def build_consensus(bamfilename, len_reference, VERBOSE=0,
                    block_len_initial=100,
                    reads_per_alignment=31,
                    store_allele_counts=False):
    '''Build a consensus from premapped and divided reads'''
    if VERBOSE:
        print 'Build consensus'

    import numpy as np
    import pysam
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet.IUPAC import ambiguous_dna
    
    from hivwholeseq.mapping_utils import align_muscle
    # Three steps:
    # 1. collect reads uniformly across the fragment
    # 2. make local consensi
    # 3. join into fragmentwide consensus
    consensi_local = []
    if store_allele_counts:
        allcounts_local = []
    pos_ref = 0
    block_len = block_len_initial
    with pysam.Samfile(bamfilename, 'rb') as bamfile:

        # Initial block
        if VERBOSE >= 2:
            print 'Block n', len(consensi_local) + 1, 
        for pos_first_block in xrange(len_reference):
            bamfile.reset()
            reads = [read for read in bamfile if read.pos == pos_first_block]
            if not len(reads):
                continue

            np.random.shuffle(reads)
            reads = reads[:n_reads_per_ali]
            seqs = [SeqRecord(Seq(read.seq[:block_len], ambiguous_dna), id=read.qname)
                    for read in reads]
            ali = np.array(align_muscle(*seqs, sort=True), 'S1', ndmin=2)
            cons_local = build_local_consensus(ali, VERBOSE=VERBOSE, store_allele_counts=store_allele_counts)
            if store_allele_counts:
                (cons_local, allcount_local) = cons_local
                allcounts_local.append(allcount_local)
            consensi_local.append(cons_local)
            pos_ref += (block_len_initial // 2) * (1 + pos_first_block // (block_len_initial // 2))
            if VERBOSE >= 2:
                print 'pos', pos_first_block, 'block len', block_len
            break

        # Add further local consensi
        while (pos_ref < len_reference):
            bamfile.reset()
            block_len = min(block_len, len_reference - pos_ref)
            if VERBOSE >= 2:
                print 'Block n', len(consensi_local) + 1, 'pos', pos_ref, 'block len', block_len
            # Get reads that cover the whole block
            reads = [read for read in bamfile
                     if (pos_ref - 100 < read.pos <= pos_ref) and \
                        (read.pos + sum(bl for (bt, bl) in read.cigar
                                        if bt in (0, 2)) >= pos_ref + block_len)]

            # Internal coverage holes are not tolerated, but the last block
            # is allowed to be missing
            if not len(reads):
                break

            np.random.shuffle(reads)
            reads = reads[:n_reads_per_ali]
            seqs = []
            for read in reads:
                pos_reft = read.pos
                # Find start of the block in the read
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
                
                seq = SeqRecord(Seq(read.seq[pos_read_start: pos_read_end],
                                    ambiguous_dna), id=read.qname)
                seqs.append(seq)

            ali = np.array(align_muscle(*seqs, sort=True), 'S1', ndmin=2)
            cons_local = build_local_consensus(ali, VERBOSE=VERBOSE, store_allele_counts=store_allele_counts)
            if store_allele_counts:
                (cons_local, allcount_local) = cons_local
                allcounts_local.append(allcount_local)
            consensi_local.append(cons_local)
            pos_ref += block_len_initial // 2

    # Join blocks (from the left)
    if VERBOSE >= 2:
        print 'Join local consensi'
    consensus = [consensi_local[0]]
    if store_allele_counts:
        allcounts = [allcounts_local[0]]
    for j, cons in enumerate(consensi_local[1:], 1):
        seed = consensus[-1][-20:]
        sl = len(seed)
        pos_start = cons.find(seed)
        # Allow imperfect matches
        if pos_start == -1:
            consm = np.fromstring(cons, 'S1')
            seedm = np.fromstring(seed, 'S1')
            n_matches = [(consm[i: i + sl] == seedm).sum()
                         for i in xrange(len(cons) + 1 - len(seed))]
            pos_start = np.argmax(n_matches)

            # Try to only add non-bogus stuff
            if n_matches[pos_start] < 0.66 * sl:
                pos_start = -1
                if VERBOSE >= 4:
                    print 'block n.', j, 'not joint!'

        if pos_start != -1:
            consensus.append(cons[pos_start + sl:])
            if store_allele_counts:
                allcounts.append(allcounts_local[j][:, pos_start + sl:])
        
    consensus = ''.join(consensus)

    if store_allele_counts:
        allcounts = np.concatenate(allcounts, axis=1)
        return (consensus, allcounts)

    return consensus




# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Build consensus, iteratively')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6i)')
    parser.add_argument('--block-length', type=int, default=100, dest='block_len',
                        help='Length of each local consensus block')
    parser.add_argument('--reads-per-alignment', type=int, default=31,
                        dest='reads_per_alignment',
                        help='Number of (random) reads used for the local consensi')
    parser.add_argument('--verbose', type=int, default=0,
                        help=('Verbosity level [0-3]'))
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--no-summary', action='store_false', dest='summary',
                        help='Do not save results in a summary file')
    parser.add_argument('--allele-counts', action='store_true', dest='allele_counts',
                        help='Also create rough allele frequencies')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    n_reads_per_ali = args.reads_per_alignment
    VERBOSE = args.verbose
    submit = args.submit
    summary = args.summary
    block_len_initial = args.block_len
    store_allele_counts = args.allele_counts

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = dataset['adapters']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # Iterate over all requested samples
    for adaID in adaIDs:

        samplename = dataset['samples'][dataset['adapters'].index(adaID)]
        # If the script is called with no fragment, iterate over all
        if not fragments:
            fragments_sample = samples[samplename]['fragments']
        else:
            from re import findall
            fragments_all = samples[samplename]['fragments']
            fragments_sample = []
            for fragment in fragments:
                frs = filter(lambda x: fragment in x, fragments_all)
                if len(frs):
                    fragments_sample.append(frs[0])

        if VERBOSE >= 3:
            print 'adaID '+adaID+': fragments '+' '.join(fragments_sample)

        for fragment in fragments_sample:

            # Submit to the cluster self if requested
            if submit:
                fork_self(seq_run, adaID, fragment, block_len_initial, n_reads_per_ali,
                          store_allele_counts, VERBOSE=VERBOSE)
                continue

            if summary:
                sfn = get_build_consensus_summary_filename(data_folder, adaID,
                                                           fragment, iterative=False)
                with open(sfn, 'w') as f:
                    f.write('Call: python build_consensus.py'+\
                            ' --run '+seq_run+\
                            ' --adaIDs '+adaID+\
                            ' --fragments '+fragment+\
                            ' --block-length '+str(block_len_initial)+\
                            ' --reads-per-alignment '+str(n_reads_per_ali)+\
                            ' --verbose '+str(VERBOSE))
                    if store_allele_counts:
                        f.write(' --allele-counts')
                    f.write('\n')

            if VERBOSE:
                print seq_run, adaID, fragment
            refseq = SeqIO.read(get_reference_premap_filename(data_folder, adaID, fragment), 'fasta')
            bamfilename = get_divided_filename(data_folder, adaID, fragment, type='bam')
            consensus = build_consensus(bamfilename, len(refseq), VERBOSE=VERBOSE,
                                        block_len_initial=block_len_initial,
                                        reads_per_alignment=n_reads_per_ali,
                                        store_allele_counts=store_allele_counts)
            if store_allele_counts:
                (consensus, allele_counts) = consensus

            # Store to file
            if VERBOSE:
                print 'Store to file'
            frag_out = fragment[:2]
            name = samplename+'_seqrun_'+seq_run+'_adaID_'+adaID+'_'+frag_out+'_consensus'
            consensusseq = SeqRecord(Seq(consensus, ambiguous_dna), id=name, name=name)

            outfile = get_consensus_filename(data_folder, adaID, frag_out, trim_primers=True)
            SeqIO.write(consensusseq, outfile, 'fasta')

            # Align all consensi via muscle and store
            ali = align_muscle(refseq, consensusseq, sort=True)
            AlignIO.write(ali, get_reference_all_filename(data_folder, adaID, fragment), 'fasta')

            if store_allele_counts:
                allele_counts.dump(get_allele_counts_filename(data_folder, adaID, frag_out))


