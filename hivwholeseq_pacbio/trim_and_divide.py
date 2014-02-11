#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/10/13
content:    Trim premapped reads according to:
            - reading into adapters (short inserts only)
            - outer and inner PCR
            - quality
            Do not trim yet for good CIGARs, i.e. indels at read edges. They
            could be simply due to distance from the reference (we are going to
            remap anyway).
'''
# Modules
import os
import argparse
from operator import itemgetter
from itertools import izip
import numpy as np
from Bio import SeqIO
import pysam

from hivwholeseq_pacbio.datasets import PacBio_runs
from hivwholeseq_pacbio.pacbio_rs_II import alpha
from hivwholeseq_pacbio.filenames import get_premapped_file, get_reference_premap_filename, \
        get_divided_filenames, get_divide_summary_filename
from hivwholeseq_pacbio.samples import samples
from hivwholeseq_pacbio.mapping_utils import test_read_integrity as test_integrity
from hivwholeseq.primer_info import primers_inner, primers_outer



# Functions
def make_output_folders(data_folder, samplename, VERBOSE=0):
    '''Make output folders'''
    from hivwholeseq.generic_utils import mkdirs
    output_filename = get_divided_filenames(data_folder, samplename, fragments=['F1'])[0]
    dirname = os.path.dirname(output_filename)
    mkdirs(dirname)
    if VERBOSE:
        print 'Folder created:', dirname


def store_reference_fragmented(data_folder, samplename, refseq, fragment_trim_poss_dict):
    '''Store FASTA files for the reference in fragments'''
    for fragment, poss in fragment_trim_poss_dict.iteritems():
        refseq_frag = refseq[poss[0]: poss[1]]
        refseq_frag.id = refseq_frag.id+'_'+fragment
        refseq_frag.name = refseq_frag.name+'_'+fragment
        refseq_frag.description = refseq_frag.description+', fragment '+fragment

        SeqIO.write(refseq_frag,
                    get_reference_premap_filename(data_folder, samplename, fragment),
                    'fasta')


def write_fragment_positions(data_folder, samplename, fragments, frags_pos):
    '''Write the fragments position to file'''
    from hivwholeseq_pacbio.filenames import get_fragment_positions_filename as gfpf
    with open(gfpf(data_folder, samplename), 'w') as f:
        f.write('\t'.join(['# Fragment', 'fwd start', 'fwd end', 'rev start', 'rev end']))
        f.write('\n')
        for (fragment, poss_full, poss_trim) in izip(fragments,
                                                     frags_pos['full'],
                                                     frags_pos['trim']):
            f.write('\t'.join(map(str, [fragment,
                                        poss_full[0], poss_trim[0],
                                        poss_trim[1], poss_full[1]]))+'\n')


def get_primer_positions(smat, fragments, type='both'):
    '''Get the primer positions for fwd, rev, or both primers'''
    from hivwholeseq.primer_info import primers_PCR
    from hivwholeseq.sequence_utils import expand_ambiguous_seq as eas

    # j controls the direction: j = 0 --> FWD, j = 1 --> REV
    types = {'fwd': [0], 'rev': [1], 'both': [0, 1]}
    js = types[type]
    primer_poss = [[], []]
    for j in js:
        pr_old_pos = 0
        pr_old = ''
        for ifr, fragment in enumerate(fragments):
            # Expand ambiguous primers in a list of all possible unambiguous ones,
            # and look for the best approximate match between all primers and a
            # sliding window in the HIV genome
            pr = primers_PCR[fragment][j]
            pr_mat = np.array(map(list, eas(pr)), 'S1')
            n_matches = [(smat[i: i + len(pr)] == pr_mat).sum(axis=1).max()
                         for i in xrange(pr_old_pos + len(pr_old),
                                         len(smat) - len(pr) + 1)]

            # NOTE: F6 rev lies in the LTR, so we risk reading the other LTR.
            # Treat it as a special case: come from the right!
            if j and ('F6' in fragment):
                pr_pos = len(smat) - len(pr) - np.argmax(n_matches[::-1])
            else:
                pr_pos = pr_old_pos + len(pr_old) + np.argmax(n_matches)

            primer_poss[j].append([pr_pos, pr_pos + len(pr)])
            pr_old_pos = pr_pos
            pr_old = pr

    if type != 'both':
        return np.array(primer_poss[type == 'rev'], int)
    else:
        return np.array(primer_poss, int)


def get_fragment_positions(smat, fragments):
    '''Get the fragment positions in the reference by approximate string match
    
    NOTE: this version of the function cannot guarantee that the end of a fragment
    comes after the start, but exploits that fragments are sorted. The bright
    side is that the function can be used to get only FWD or REV primers, ignoring
    the rest of the output.
    '''
    primer_poss = get_primer_positions(smat, fragments, type='both')

    fragment_full_poss = [primer_poss[0, :, 0],
                          primer_poss[1, :, 1]]
    fragment_trim_poss = [primer_poss[0, :, 1],
                          primer_poss[1, :, 0]]

    return {'full': np.array(fragment_full_poss, int).T,
            'trim': np.array(fragment_trim_poss, int).T}


def test_fragment_assignment(read, n_frag, len_frag, VERBOSE=True):
    '''Test assignment of the read pair to a fragment'''

    if (read.pos <  0) or (read.pos + abs(read.isize) > len_frag):
        if VERBOSE:
            print 'Read pair is misassigned to fragment:', read.qname, n_frag
        return True

    return False


def test_sanity(read, n_frag, len_frag):
    '''Perform sanity checks on supposedly good reads'''
    if not test_integrity(read):
        return False

    if not test_fragment_assignment(read, n_frag, len_frag):
        return False

    return True


def test_outer_primer(read, primers_out_pos, primers_out_seq, lref):
    '''Test for reads starting at outer primer.'''
    # The first fwd and last rev outer primers are special: we can spot them for
    # it is mapped at 0 with an insert as first CIGAR (a read never starts with
    # a deletion) or it reaches the end of the reference with a final insert.
    if (read.pos == 0) and (read.cigar[0][0] == 1):
        return True
    if (read.pos + abs(read.isize) == lref) and (read.cigar[-1][0] == 1):
        return True

    # TODO: do you work in case there is only fwd OR rev?
    # The internal outer primer cannot be checked by positions only, because
    # there might be little indels in the reads and we want to be conservative;
    # however, the read start/end must be mapped somewhere close
    # FWD primer
    for pr_pos, pr_mat in izip(primers_out_pos['fwd'], primers_out_seq['fwd']):
        if np.abs(read.pos - pr_pos) < 8:
            rfwd = np.fromstring(read.seq[:pr_mat.shape[1]], 'S1')
            if (rfwd == pr_mat).mean(axis=1).max() > 0.9:
                return True
            else:
                # If the read is close to an outer primer, it is far from any other
                break

    # REV primer
    for primer_pos, pr_mat in izip(primers_out_pos['rev'], primers_out_seq['rev']):
        if np.abs(read.pos + abs(read.isize) - pr_pos) < 8:
            rrev = np.fromstring(read.seq[-pr_mat.shape[1]:], 'S1')
            if (rrev == pr_mat).mean(axis=1).max() > 0.9:
                return True
            else:
                break

    return False


def assign_to_fragment(read, fragment_full_poss):
    '''Assign read pair to fragments'''

    # Insert coordinates
    ins_start = read.pos
    ins_end = read.pos + abs(read.isize)

    # Fragment coordinates (including primers)
    fragment_start, fragment_end = fragment_full_poss.T

    # Check the insert
    frags_ins = ((ins_start >= fragment_start) &
                 (ins_start < fragment_end) &
                 (ins_end <= fragment_end) \
                ).nonzero()[0]

    return frags_ins


def trim_primers(read, frag_pos, include_tests=False):
    '''Trim inner primers
    
    - frag_pos are the coordinates of the fragment trimmed of innermost primers
    '''

    if include_tests:
        if test_integrity(read):
            print 'trim_primers (entry):'
            import ipdb; ipdb.set_trace()

    # FWD primer
    if read.pos < frag_pos[0]:
        ref_pos = read.pos
        read_pos = 0
        cigar = list(read.cigar)
        for (bt, bl) in read.cigar:
            if bt == 0:
                # Strictly greater, because we dump the CIGAR ending at the
                # right position. In the corner case, ref_pos == frag_pos[0]
                # and we add the whole block (and do not move read_pos).
                if ref_pos + bl > frag_pos[0]:
                    cigar[0] = (bt, ref_pos + bl - frag_pos[0])
                    read_pos += frag_pos[0] - ref_pos
                    ref_pos = frag_pos[0]
                    break
                else:
                    cigar.pop(0)
                    read_pos += bl
                    ref_pos += bl

            # Forget insertions right after the primer
            elif bt == 1:
                cigar.pop(0)
                read_pos += bl

            # Starting with a deletion is not allowed
            elif bt == 2:
                cigar.pop(0)
                ref_pos += bl
                if ref_pos > frag_pos[0]:
                    break

        # if we chopped off everything, trash
        if not len(cigar):
            return True

        seq = read.seq
        qual = read.qual
        read.pos = ref_pos
        read.seq = seq[read_pos:]
        read.qual = qual[read_pos:]
        read.cigar = cigar
        isize = sum(bl for (bt, bl) in read.cigar if bt in (0, 2))
        if read.is_reverse:
            read.isize = -isize
        else:
            read.isize = isize

    # REV primer (come from the right)
    ref_pos = read.pos + abs(read.isize)
    if ref_pos > frag_pos[1]:
        read_pos = read.rlen
        cigar = list(read.cigar[::-1])
        for (bt, bl) in read.cigar[::-1]:
            if bt == 0:
                # Strictly less, because we dump the CIGAR starting at the
                # right position. In the corner case, ref_pos == frag_pos[1]
                # and we add the whole block (and do not move read_pos).
                if ref_pos - bl < frag_pos[1]:
                    cigar[0] = (bt, frag_pos[1] - (ref_pos - bl))
                    read_pos -= ref_pos - frag_pos[1]
                    break
                else:
                    cigar.pop(0)
                    read_pos -= bl
                    ref_pos -= bl
            
            # Forget insertions right before the rev primer
            elif bt == 1:
                cigar.pop(0)
                read_pos -= bl
            
            # Ending with a deletion is not allowed
            elif bt == 2:
                cigar.pop(0)
                ref_pos -= bl
                if ref_pos < frag_pos[1]:
                    break

        # if we chopped off everything, trash
        if not len(cigar):
            return True

        seq = read.seq
        qual = read.qual
        read.seq = seq[:read_pos]
        read.qual = qual[:read_pos]
        read.cigar = cigar[::-1]
        isize = sum(bl for (bt, bl) in read.cigar if bt in (0, 2))
        if read.is_reverse:
            read.isize = -isize
        else:
            read.isize = isize

    if include_tests:
        if test_integrity(read):
            print 'trim_primers (exit):'
            import ipdb; ipdb.set_trace()

    return False


def trim_and_divide_reads(data_folder, samplename, fragments,
                          maxreads=-1, VERBOSE=0,
                          minisize=100,
                          include_tests=False, summary=True):
    '''Trim reads and divide them into fragments'''
    if VERBOSE:
        print 'Trim and divide into fragments: sample name '+samplename+', fragments: '+\
                ' '.join(fragments)

    if summary:
        with open(get_divide_summary_filename(data_folder, samplename), 'a') as f:
            f.write('Fragments used: '+' '.join(fragments)+'\n')

    ref_filename = get_reference_premap_filename(data_folder, samplename)
    refseq = SeqIO.read(ref_filename, 'fasta')
    smat = np.array(refseq, 'S1')
    len_reference = len(refseq)

    # Get the positions of fragment start/end, w/ and w/o primers
    frags_pos = get_fragment_positions(smat, fragments)
    store_reference_fragmented(data_folder, samplename, refseq,
                               dict(zip(*[fragments, frags_pos['trim']])))
    if summary:
        with open(get_divide_summary_filename(data_folder, samplename), 'a') as f:
            f.write('Primer positions (for fragments):\n')
            for (fragment, poss_full, poss_trim) in izip(fragments,
                                                         frags_pos['full'],
                                                         frags_pos['trim']):
                f.write(fragment+': fwd '+str(poss_full[0])+' '+str(poss_trim[0])+\
                                 ', rev '+str(poss_trim[1])+' '+str(poss_full[1])+'\n')
    write_fragment_positions(data_folder, samplename, fragments, frags_pos)

    # Get the positions of the unwanted outer primers (in case we DO nested PCR
    # for that fragment)
    # NOTE: the LTRs make no problem, because the rev outer primer of F6
    # is not in the reference anymore if F6 has undergone nested PCR
    from re import findall
    primers_out = {'fwd': [], 'rev': []}
    for i, fr in enumerate(fragments):
        if (i != 0) and findall(r'F[2-6][a-z]?i', fr):
            primers_out['fwd'].append(fr[:-1]+'o')
        if (i != len(fragments) - 1) and findall(r'F[1-5][a-z]?i', fr):
            primers_out['rev'].append(fr[:-1]+'o')

    # Get all possible unambiguous primers for the unwanted outer primers
    from hivwholeseq.primer_info import primers_PCR
    from hivwholeseq.sequence_utils import expand_ambiguous_seq as eas
    primers_out_seq = {'fwd': [np.array(map(list, eas(primers_PCR[fr][0])),
                                        'S1', ndmin=2)
                               for fr in primers_out['fwd']],
                       'rev': [np.array(map(list, eas(primers_PCR[fr][1])),
                                        'S1', ndmin=2)
                               for fr in primers_out['rev']],
                      }

    primers_out_pos = {'fwd': [], 'rev': []}
    if primers_out['fwd']:
        primers_out_pos['fwd'] = get_primer_positions(smat, primers_out['fwd'], 'fwd')[:, 0]
    if primers_out['rev']:
        primers_out_pos['rev'] = get_primer_positions(smat, primers_out['rev'], 'rev')[:, 1]

    # Input and output files
    input_filename = get_premapped_file(data_folder, samplename, type='bam')
    if not os.path.isfile(input_filename):
        convert_sam_to_bam(input_filename)
    output_filenames = get_divided_filenames(data_folder, samplename, fragments, type='bam')
    with pysam.Samfile(input_filename, 'rb') as bamfile:
        with pysam.Samfile(output_filenames[0], 'wb', template=bamfile) as fo_F1,\
             pysam.Samfile(output_filenames[1], 'wb', template=bamfile) as fo_F2,\
             pysam.Samfile(output_filenames[2], 'wb', template=bamfile) as fo_F3,\
             pysam.Samfile(output_filenames[3], 'wb', template=bamfile) as fo_F4,\
             pysam.Samfile(output_filenames[4], 'wb', template=bamfile) as fo_F5,\
             pysam.Samfile(output_filenames[5], 'wb', template=bamfile) as fo_F6,\
             pysam.Samfile(output_filenames[6], 'wb', template=bamfile) as fo_am,\
             pysam.Samfile(output_filenames[7], 'wb', template=bamfile) as fo_cm,\
             pysam.Samfile(output_filenames[8], 'wb', template=bamfile) as fo_um,\
             pysam.Samfile(output_filenames[9], 'wb', template=bamfile) as fo_lq:

            # Collect the file handles
            file_handles = (fo_F1, fo_F2, fo_F3, fo_F4, fo_F5, fo_F6)

            # Iterate over the mapped reads and assign fragments
            n_mapped = [0 for fragment in fragments]
            n_unmapped = 0
            n_crossfrag = 0
            n_ambiguous = 0
            n_outer = 0
            n_lowq = 0
            for irp, read in enumerate(bamfile):

                if irp == maxreads:
                    if VERBOSE:
                        print 'Maximal number of reads reached:', maxreads
                    break

                if VERBOSE >= 2:
                    if not ((irp+1) % 10000):
                        print irp+1, maxreads

                # If unmapped or mini, discard
                if read.is_unmapped or (read.rlen < 50):
                    if VERBOSE >= 3:
                        print 'Read unmapped/tiny:', read.qname
                    n_unmapped += 1
                    fo_um.write(read)
                    continue

                # If the insert is a misamplification from the outer primers
                # in fragments that underwent nested PCR,
                # trash it (it will have skewed amplification anyway). We cannot
                # find all of those, rather only the ones still carrying the
                # primer itself (some others have lost it while shearing). For
                # those, no matter what happens at the end (reading into adapters,
                # etc.), ONE of the reads in the pair will start exactly with one
                # outer primer: if the rev read with a rev primer, if the fwd
                # with a fwd one. Test all six.
                if (len(primers_out_pos['fwd']) or len(primers_out_pos['rev'])) and \
                   test_outer_primer(read,
                                     primers_out_pos, primers_out_seq,
                                     len_reference):
                    if VERBOSE >= 3:
                        print 'Read from outer primer:', read.qname
                    n_outer += 1
                    fo_um.write(read)
                    continue

                # Assign to a fragment now, so that primer trimming is faster 
                frags_read = assign_to_fragment(read, frags_pos['full'])

                # 1. If no fragments are possible (e.g. one read crosses the
                # fragment boundary, they map to different fragments), dump it
                # into a special bucket
                if len(frags_read) == 0:
                    n_crossfrag += 1
                    fo_cm.write(read)
                    continue

                # 2. If 2+ fragments are possible (tie), put into a special bucket
                # (essentially excluded, because we want two independent measurements
                # in the overlapping region, but we might want to recover them)
                elif len(frags_read) > 1:
                    n_ambiguous += 1
                    fo_am.write(read)
                    continue

                # 3. If the intersection is a single fragment, good: trim the primers
                n_frag = frags_read[0]
                trashed_primers = trim_primers(read, frags_pos['trim'][n_frag],
                                               include_tests=include_tests)
                if trashed_primers:
                    if VERBOSE >= 3:
                        print 'Read mismapped:', read.qname
                    n_unmapped += 1
                    fo_um.write(read)
                    continue

                # TODO: Quality trimming: not implemented yet

                # Change coordinates into the fragmented reference (primer-trimmed)
                read.pos -= frags_pos['trim'][n_frag, 0]

                if include_tests:
                    lfr = frags_pos['trim'][n_frag, 1] - frags_pos['trim'][n_frag, 0]
                    if test_sanity(read, n_frag, lfr):
                        print 'Tests failed:', read.qname
                        import ipdb; ipdb.set_trace()

                # Store to file the good reads
                n_frag_all = int(fragments[n_frag][1]) - 1
                n_mapped[n_frag_all] += 1
                file_handles[n_frag_all].write(read)

    if VERBOSE:
        print 'Trim and divide results: sample name '+samplename
        print 'Total:\t\t', irp
        print 'Mapped:\t\t', sum(n_mapped), n_mapped
        print 'Unmapped:\t', n_unmapped
        print 'Outer primer\t', n_outer
        print 'Crossfrag:\t', n_crossfrag
        print 'Ambiguous:\t', n_ambiguous
        print 'Low-quality:\t', n_lowq

    # Write summary to file
    if summary:
        with open(get_divide_summary_filename(data_folder, samplename), 'a') as f:
            f.write('\n')
            f.write('Trim and divide results: sample name '+samplename+'\n')
            f.write('Total:\t\t'+str(irp)+'\n')
            f.write('Mapped:\t\t'+str(sum(n_mapped))+' '+str(n_mapped)+'\n')
            f.write('Unmapped:\t'+str(n_unmapped)+'\n')
            f.write('Outer primer\t'+str(n_outer)+'\n')
            f.write('Crossfrag:\t'+str(n_crossfrag)+'\n')
            f.write('Ambiguous:\t'+str(n_ambiguous)+'\n')
            f.write('Low-quality:\t'+str(n_lowq)+'\n')



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Trim and divide reads into fragments')
    parser.add_argument('--run', default='Upp23',
                        help='PacBio run to analyze (e.g. Upp23)')
    parser.add_argument('--samples', required=True, nargs='*',
                        help='Samples to analyze (e.g. S1 S2)')
    parser.add_argument('--verbose', default=0, type=int,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to analyze')
    parser.add_argument('--minisize', type=int, default=100,
                        help='Minimal insert size to keep')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--test', action='store_true',
                        help='Include sanity checks on mapped reads (slow)')
    parser.add_argument('--no-summary', action='store_false', dest='summary',
                        help='Do not save results in a summary file')

    args = parser.parse_args()
    seq_run = args.run
    samplenames = args.samples
    VERBOSE = args.verbose
    maxreads = args.maxreads
    minisize = args.minisize
    submit = args.submit
    include_tests = args.test
    summary = args.summary

    # Specify the dataset
    dataset = PacBio_runs[seq_run]
    data_folder = dataset['folder']

    if not samplenames:
        samplenames = PacBio_runs[seq_run]['samples']

    for samplename in samplenames:

        # Submit to the cluster self if requested
        if submit:
            if include_tests:
                raise ValueError('Tests require an interactive shell')
            fork_self(seq_run, samplename, VERBOSE=VERBOSE, maxreads=maxreads,
                      minisize=minisize, summary=summary)
            continue

        make_output_folders(data_folder, samplename, VERBOSE=VERBOSE)

        if summary:
            with open(get_divide_summary_filename(data_folder, samplename), 'w') as f:
                f.write('Call: python trim_and_divide.py --run '+seq_run+\
                        ' --samples '+samplename+\
                        ' --minisize '+str(minisize)+\
                        ' --verbose '+str(VERBOSE))
                if maxreads != -1:
                    f.write(' --maxreads '+str(maxreads))
                if include_tests:
                    f.write(' --include_tests')
                f.write('\n')

        fragments = samples[samplename]['fragments']

        # Trim reads and assign them to a fragment (or discard)
        # Note: we pass over the file only once
        trim_and_divide_reads(data_folder, samplename, fragments,
                              maxreads=maxreads, VERBOSE=VERBOSE,
                              minisize=minisize,
                              include_tests=include_tests,
                              summary=summary)
        continue
