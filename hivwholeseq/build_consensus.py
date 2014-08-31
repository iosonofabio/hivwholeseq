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

from hivwholeseq.samples import load_sequencing_run, SampleSeq
from hivwholeseq.mapping_utils import align_muscle
from hivwholeseq.filenames import get_divided_filename, \
        get_premapped_filename, \
        get_reference_premap_filename, \
        get_consensus_filename, \
        get_allele_counts_filename, \
        get_build_consensus_summary_filename, \
        get_reference_consensus_ali_filename
from hivwholeseq.fork_cluster import fork_build_consensus as fork_self



# Functions
def build_local_consensus(seqs, VERBOSE=0, store_allele_counts=False, full_cover=True):
    '''Build a local consensus from an MSA'''
    # There is only ONE tricky point: what to do if some reads do not cover the whole
    # block, e.g. at the end of a fragment because of low coverage?
    # If full_cover == False, convert MSA gaps at the end of too short reads into N

    import numpy as np
    from hivwholeseq.miseq import alpha
    from hivwholeseq.mapping_utils import align_muscle

    ali = np.array(align_muscle(*seqs, sort=True), 'S1', ndmin=2)
    if full_cover:
        allele_counts = np.array([(ali == a).sum(axis=0) for a in alpha], int, ndmin=2)
    else:
        allele_counts = np.zeros((len(alpha), len(ali[0])),int)
        for i in xrange(len(seqs)):
            if ali[i, -1] == '-':
                first_finalgap = len(ali[i].tostring().rstrip('-'))
                ali[i, first_finalgap:] = 'X'
            for ai, a in enumerate(alpha):
                allele_counts[ai] += ali[i] == a

        cov = allele_counts.sum(axis=0)
        allele_counts = allele_counts[:, cov > 0]

    cons_local = []
    for counts in allele_counts.T:
        # Pick max count nucleotide, ignoring N
        maxinds = (counts[:-1] == counts.max()).nonzero()[0]
        if len(maxinds) < 1:
            cons_local.append('-')
            continue
        # Pick a random nucleotide in case of a tie
        elif len(maxinds) > 1:
            np.random.shuffle(maxinds)
        maxind = maxinds[0]
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
                    accept_holes=False,
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
    consensus = None
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

            # The first block has to make a consensus for the FIRST base, this needs
            # at least ONE read starting exactly at the first position. Otherwise,
            # the same is repeated for position 2, and so on.
            reads = [read for read in bamfile if (read.is_proper_pair) and (read.pos == pos_first_block)]
            if not len(reads):
                continue

            np.random.shuffle(reads)
            reads = reads[:n_reads_per_ali]
            seqs = [SeqRecord(Seq(read.seq[:block_len], ambiguous_dna), id=read.qname)
                    for read in reads]
            cons_local = build_local_consensus(seqs, VERBOSE=VERBOSE, store_allele_counts=store_allele_counts)
            if store_allele_counts:
                (cons_local, allcount_local) = cons_local
                allcounts_local.append(allcount_local)
            consensi_local.append(cons_local)
            pos_ref += (block_len_initial // 2) * (1 + pos_first_block // (block_len_initial // 2))
            if VERBOSE >= 2:
                print 'pos', pos_first_block, 'to', pos_first_block + block_len, 'block len', block_len
            break

        # Start consensus
        if len(consensi_local) == 1:
            consensus = [consensi_local[0]]
            if store_allele_counts:
                allcounts = [allcounts_local[0]]

        # Divide reads by block (more efficient than scrolling the file every time)
        # FIXME: extract random subsample, assign to blocks, and only complete the missing blocks!
        reads_by_block = [[] for n_block in xrange((len_reference - pos_ref) // (block_len_initial // 2))]
        bamfile.reset()
        for read in bamfile:
            if not read.is_proper_pair:
                continue
            pos_ref_tmp = pos_ref
            n_block = 1
            while (pos_ref_tmp < len_reference):
                block_len_tmp = min(block_len, len_reference - pos_ref)
                read_start = read.pos
                read_end = read.pos + sum(bl for (bt, bl) in read.cigar if bt in (0, 2))
                if (pos_ref_tmp - 100 < read_start <= pos_ref_tmp) and \
                   (read_end >= pos_ref_tmp + block_len_tmp):
                    reads_by_block[n_block - 1].append(read)
                    break

                pos_ref_tmp += block_len_initial // 2
                n_block += 1

        # Stack local consensi on top of the first one
        n_block = 1
        while (pos_ref < len_reference):
            block_len = min(block_len, len_reference - pos_ref)
            if block_len < block_len_initial // 2:
                break
            if VERBOSE >= 2:
                print 'Block n', len(consensi_local) + 1, 'pos', pos_ref, 'to', pos_ref + block_len, 'block len', block_len

            # Get reads that cover the whole block
            reads = reads_by_block[n_block - 1]
            n_block += 1

            #FIXME
            #if n_block >= 2:
            #    print pos_ref, pos_ref + block_len
            #    import ipdb; ipdb.set_trace()

            # Internal coverage holes are not tolerated, but the last block
            # is allowed to be missing. However, we should try to squeeze out
            # all the bases by rescanning the reads a last time with less strict
            # criteria: if it has even one base more than what we have, add it
            if len(reads):
                full_cover= True
            else:
                full_cover= False
                bamfile.reset()
                reads = []
                for read in bamfile:
                    if not read.is_proper_pair:
                        continue
                    read_start = read.pos
                    read_end = read.pos + sum(bl for (bt, bl) in read.cigar if bt in (0, 2))
                    if (read_start <= pos_ref) and (read_end > pos_ref + block_len_initial // 2):
                        reads.append(read)

                if not len(reads):
                    if pos_ref + block_len < len_reference:
                        if VERBOSE >= 2:
                            print 'WARNING: consensus looks interrupted in mid-way'
                    break

            # Take a random subsample of reads. If it's a problematic block, not
            # fully covered, take more reads than usual
            if full_cover:
                np.random.shuffle(reads)
                reads = reads[:n_reads_per_ali]
            else:
                # Trim all, then take longest
                pass

            # Trim reads from the left to start all at the block start
            # NOTE: reads have been selected to start @ or before the block start!
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

            # If it's a problematic block, take longest reads
            if not full_cover:
                seqs.sort(key=len, reverse=True)
                seqs = seqs[:n_reads_per_ali]

            #FIXME
            #if n_block >= 2:
            #    print pos_ref, pos_ref + block_len
            #    import ipdb; ipdb.set_trace()

            # Make local consensus
            cons_local = build_local_consensus(seqs, VERBOSE=VERBOSE,
                                               store_allele_counts=store_allele_counts,
                                               full_cover=full_cover)
            if store_allele_counts:
                (cons_local, allcount_local) = cons_local
                allcounts_local.append(allcount_local)
            consensi_local.append(cons_local)

            pos_ref += block_len_initial // 2

            # Join block <-- to the stack
            if consensus is None:
                consensus = [consensi_local[0]]
                if store_allele_counts:
                    allcounts = [allcounts_local[0]]
            else:
                cons = cons_local
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
                            print 'Block n.', len(consensi_local)+': cannot stack to previous one!'
        
                if pos_start != -1:
                    consensus.append(cons[pos_start + sl:])
                    if store_allele_counts:
                        allcounts.append(allcounts_local[-1][:, pos_start + sl:])
                
                elif accept_holes:
                    consensus.append('N' * 10)
                    consensus.append(cons)
                    if store_allele_counts:
                        tmpall = np.zeros((allcounts_local[-1].shape[0], 10), int)
                        tmpall[-1] = 1
                        allcounts.append(tmpall)
                        allcounts.append(allcounts_local[-1])

    if consensus is None:
        raise ValueError('Consensus is still None: unable to build!')

    consensus = ''.join(consensus)

    if store_allele_counts:
        allcounts = np.concatenate(allcounts, axis=1)
        return (consensus, allcounts)

    return consensus




# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Build consensus by mapping-assisted assembly',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6 genomewide)')
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
    dataset = load_sequencing_run(seq_run)
    data_folder = dataset.folder

    # If the script is called with no adaID, iterate over all
    dataset.discard_nondivided_samples()
    samples = dataset.samples
    if adaIDs is not None:
        samples = samples.loc[samples.adapter.isin(adaIDs)]

    if VERBOSE >= 3:
        print 'adaIDs', samples.adapter

    # Iterate over all requested samples
    for samplename, sample in samples.iterrows():
        sample = SampleSeq(sample)
        adaID = sample.adapter

        # If the script is called with no fragment, iterate over all
        if not fragments:
            fragments_sample = sample.regions_complete
        else:
            from re import findall
            fragments_all = sample.regions_complete
            fragments_sample = []
            for fragment in fragments:
                frs = filter(lambda x: fragment in x, fragments_all)
                if len(frs):
                    fragments_sample.append(frs[0])
            if 'genomewide' in fragments:
                fragments_sample.append('genomewide')

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
            if fragment == 'genomewide':
                refseq = SeqIO.read(get_reference_premap_filename(data_folder, adaID), 'fasta')
                bamfilename = get_premapped_filename(data_folder, adaID, type='bam')
                frag_out = fragment
            else:
                fn = get_reference_premap_filename(data_folder, adaID, fragment)
                bamfilename = get_divided_filename(data_folder, adaID, fragment, type='bam')

                #FIXME: old nomenclature for F3a is F3
                if not os.path.isfile(fn) and fragment[:3] == 'F3a':
                    fn = get_reference_premap_filename(data_folder, adaID, 'F3'+fragment[-1])
                if not os.path.isfile(bamfilename) and fragment[:3] == 'F3a':
                    bamfilename = get_divided_filename(data_folder, adaID, 'F3'+fragment[-1], type='bam')

                refseq = SeqIO.read(fn, 'fasta')
                frag_out = fragment[:2]

            consensus = build_consensus(bamfilename, len(refseq), VERBOSE=VERBOSE,
                                        block_len_initial=block_len_initial,
                                        reads_per_alignment=n_reads_per_ali,
                                        accept_holes=(fragment == 'genomewide'),
                                        store_allele_counts=store_allele_counts)
            if store_allele_counts:
                (consensus, allele_counts) = consensus

            # Store to file
            if VERBOSE:
                print 'Store to file'
            name = samplename+'_seqrun_'+seq_run+'_adaID_'+adaID+'_'+frag_out+'_consensus'
            consensusseq = SeqRecord(Seq(consensus, ambiguous_dna), id=name, name=name)

            # Align consensus to reference via muscle and trim end gaps in ref
            # (improper primer trimming in trim_and_divide)
            ali = align_muscle(refseq, consensusseq, sort=True)

            if ali[0][-1] == '-':
                start_nongap = len(ali[0]) - len(ali[0].seq.lstrip('-'))
                end_nongap = len(ali[0].seq.rstrip('-'))
                ali = ali[:, start_nongap: end_nongap]

            if VERBOSE >= 2:
                print ali[:, :30]
                print ali[:, -30:]
                print 'Lenghts: ref', len(refseq), 'consensus', len(consensusseq)
                len_ali = ali.get_alignment_length()
                n_diff = sum(ali[0, i] != ali[1, i] for i in xrange(len_ali))
                print 'Differences from ref:', n_diff, '('+'{:3.1f}'.format(100.0 * n_diff / len_ali)+'%)'

            # Ungap consensus
            consensusseq = SeqRecord(ali[1].seq, id=name, name=name)
            if '-' in consensusseq:
                consensusseq.seq = consensusseq.seq.ungap('-')

            # Write output
            outfile = get_consensus_filename(data_folder, adaID, frag_out, trim_primers=True)
            SeqIO.write(consensusseq, outfile, 'fasta')

            AlignIO.write(ali, get_reference_consensus_ali_filename(data_folder, adaID, fragment), 'fasta')

            if store_allele_counts:
                allele_counts.dump(get_allele_counts_filename(data_folder, adaID, frag_out))


