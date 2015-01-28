# vim: fdm=marker
'''
author:     Fabio Zanini
date:       13/09/14
content:    Examine reads that are more than the threshold from consensus
            (typically 30 changes) to find out what they are.
'''
# Modules
import sys
import os
import argparse
from operator import itemgetter
from itertools import izip
import numpy as np
import pysam
from seqanpy import align_global, align_local, align_overlap
from Bio import AlignIO

from hivwholeseq.sequencing.filenames import get_filter_mapped_summary_filename, \
        get_mapped_filename
from hivwholeseq.sequencing.samples import SampleSeq, load_samples_sequenced
from hivwholeseq.patients.filenames import get_consensi_alignment_filename
from hivwholeseq.utils.generic import getchar
from hivwholeseq.utils.sequence import pretty_print_pairwise_ali



# Functions
def fish_distant_reads(bamfilename, ref,
                       min_mismatches=20, max_mismatches=30,
                       VERBOSE=0, maxseqs=-1):
    '''Fish distant reads from the trash'''
    import numpy as np

    from hivwholeseq.utils.mapping import pair_generator, reads_to_seqrecord
    from hivwholeseq.sequencing.filter_mapped_reads import check_overhanging_reads, \
            get_distance_from_consensus

    distances = []
    seqs = []
    edges = []
    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        for irp, reads in enumerate(pair_generator(bamfile)):
            if VERBOSE >= 2:
                if not ((irp + 1) % 10000):
                    print irp + 1

            (read1, read2) = reads
            i_fwd = reads[0].is_reverse

            # Check a few things to make sure we are looking at paired reads
            if read1.qname != read2.qname:
                raise ValueError('Read pair '+str(irp)+': reads have different names!')

            # Ignore unmapped reads
            if read1.is_unmapped or read2.is_unmapped:
                continue
            
            # Ignore not properly paired reads (this includes mates sitting on
            # different fragments)
            if (not read1.is_proper_pair) or (not read2.is_proper_pair):
                continue

            # Check for overhangs beyond the edge
            skip = check_overhanging_reads(reads, len(ref))
            if skip:
                continue

            # Fish out our reads
            dc = get_distance_from_consensus(ref, reads, VERBOSE=VERBOSE)
            if (min_mismatches <= dc.sum() <= max_mismatches):
                if VERBOSE >= 3:
                    print 'Gotcha!', reads[0].qname
                seqs.append(reads[0])
                seqs.append(reads[1])
                distances.append(dc)
                edge = [(read.pos, read.pos + sum(bl for bt, bl in read.cigar if bt in (0, 2)))
                        for read in reads]
                edges.append(edge)

                if len(seqs) // 2 == maxseqs:
                    if VERBOSE >= 2:
                        print 'Max seqs reached:', maxseqs
                    break

        seqs = list(pair_generator(reads_to_seqrecord(seqs)))

    distances = np.array(distances, int)
    return (distances, edges, seqs)



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Examine distant reads',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    runs_or_samples = parser.add_mutually_exclusive_group(required=True)
    runs_or_samples.add_argument('--runs', nargs='+',
                                 help='Seq runs to analyze (e.g. Tue28, test_tiny)')
    runs_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to map (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--adaIDs', nargs='+',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--fragments', nargs='+',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--nopatients', action='store_false', dest='use_pats',
                        help='Include non-patient samples (e.g. reference strains)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxseqs', type=int, default=-1,
                        help='Stop when collected so many seqs')
    parser.add_argument('--min-mismatches', type=int, default=20,
                        dest='min_mismatches',
                        help='Minimal number of mismatches to select')
    parser.add_argument('--max-mismatches', type=int, default=20,
                        dest='max_mismatches',
                        help='Maximal number of mismatches to select')

    args = parser.parse_args()
    seq_runs = args.runs
    samplenames = args.samples
    adaIDs = args.adaIDs
    fragments = args.fragments
    use_pats = args.use_pats
    VERBOSE = args.verbose
    maxseqs = args.maxseqs
    min_mismatches = args.min_mismatches
    max_mismatches = args.max_mismatches

    samples = load_samples_sequenced()
    if seq_runs is not None:
        samples = samples.loc[samples['seq run'].isin(seq_runs)]
    
        if adaIDs is not None:
            samples = samples.loc[samples.adapter.isin(adaIDs)]
    
        if use_pats:
            samples = samples.loc[samples['patient sample'] != 'nan']
    else:
        samples = samples.loc[samplenames]

    if fragments is None:
        fragments = ['F'+str(i+1) for i in xrange(6)]

    alis = {fr: AlignIO.read(get_consensi_alignment_filename('all', fr), 'fasta')
            for fr in fragments}

    for samplename, sample in samples.iterrows():
        sample = SampleSeq(sample)
        data_folder = sample.seqrun_folder
        adaID = sample.adapter
        pname = sample.patientname

        for fragment in fragments:
            if VERBOSE >= 1:
                print sample['seq run'], adaID, fragment, samplename,

            # Read the summary filename of the filter_mapped, and find out whether
            # there are many distant reads (a few are normal)
            fn = get_filter_mapped_summary_filename(data_folder, adaID, fragment)
            if os.path.isfile(fn):
                found = False
                with open(fn, 'r') as f:
                    for line in f:
                        line = line.rstrip('\n')
                        if line[:4] == 'Good':
                            n_good = int(line.split()[-1])

                        elif line[:14] == 'Many-mutations':
                            n_distant = int(line.split()[-1])
                            found = True
                            break

                if not found:
                    if VERBOSE >= 1:
                        print 'not filtered (probably no HIV reads)'
                    continue

                frac_dist = 1.0 * n_distant / n_good
                if frac_dist < 0.01:
                    if VERBOSE >= 1:
                        print '< 1% of reads are distant'
            
                else:
                    if VERBOSE >= 1:
                        print '{:3.0%}'.format(frac_dist), 'of reads are distant'

            consrec = sample.get_consensus(fragment)
            bamfilename = get_mapped_filename(data_folder, adaID, fragment,
                                              filtered=False)

            (ds, edges, seqs) = fish_distant_reads(bamfilename, consrec, VERBOSE=VERBOSE,
                                                   min_mismatches=min_mismatches,
                                                   max_mismatches=max_mismatches,
                                                   maxseqs=maxseqs)
            indrandom = np.arange(len(ds))
            np.random.shuffle(indrandom)
            ds = ds[indrandom]
            edges = np.array(edges)[indrandom]
            seqs = [seqs[i] for i in indrandom]

            for irp, (dpair, edgepair, seqpair) in enumerate(izip(ds, edges, seqs)):
                # NOTE: Take only the most distant read of a pair
                print irp, dpair

                i = dpair.argmax()
                d = dpair[i]
                edge = edgepair[i]
                seq = seqpair[i]

                (score, ali1, ali2) = align_global(seq, consrec[edge[0]: edge[1]])
                scoremax = 3 * len(ali1)
                delta = scoremax - score
                ali = [ali2, ali1]

                print 'Alignment to its own consensus (delta = '+str(delta)+')'
                pretty_print_pairwise_ali(ali,
                                          'cons',
                                          'read'+str(i+1)+' '+str(edge),
                                          len_name=25, width=90)
                print ''

                # Compare to all consensi and find the closest
                alifr = alis[fragment]
                alifrpw = []
                for cons in alifr:
                    alifrpw.append(align_overlap(cons.seq.ungap('-'), seq))
                scores = map(itemgetter(0), alifrpw)
                indmax = np.argmax(scores)

                alimax = alifrpw[indmax][1:]
                start = len(alimax[1]) - len(alimax[1].lstrip('-'))
                end = len(alimax[1].rstrip('-'))
                alimax = [s[start: end] for s in alimax]
                score = scores[indmax]
                scoremax = 3 * len(alimax[0])
                delta = scoremax - score

                name1 = ' '.join(['cons '] + alifr[indmax].name.split('_')[::2])
                name2 = ' '.join(['read'+str(i+1), pname, sample['patient sample']])
                print 'Alignment to best consensus (delta = '+str(delta)+')'
                pretty_print_pairwise_ali(alimax, name1, name2, len_name=25, width=90)
                print ''

                print 'Press q to exit'
                ch = getchar()
                if ch.lower() == 'q':
                    break 

            break
