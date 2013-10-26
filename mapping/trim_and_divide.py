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
import subprocess as sp
import argparse
import numpy as np
import pysam

from mapping.miseq import alpha
from mapping.datasets import MiSeq_runs
from mapping.primer_info import primers_coordinates_HXB2_inner as pcis
from mapping.primer_info import primers_coordinates_HXB2_outer as pcos
from mapping.primer_info import primers_inner, primers_outer
from mapping.filenames import get_HXB2_entire, get_premapped_file, \
        get_divided_filenames, get_divide_summary_filename
from mapping.mapping_utils import pair_generator, convert_sam_to_bam


# Globals
# Cluster submit
import mapping
JOBDIR = mapping.__path__[0].rstrip('/')+'/'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
JOBSCRIPT = JOBDIR+'trim_and_divide.py'
cluster_time = '2:59:59'
vmem = '1G'



# Functions
def fork_self(miseq_run, adaID, VERBOSE=0, maxreads=-1):
    '''Fork self for each adapter ID'''
    if VERBOSE:
        print 'Forking to the cluster: adaID '+'{:02d}'.format(adaID)

    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'tr+div '+'{:02d}'.format(adaID),
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', miseq_run,
                 '--adaIDs', adaID,
                 '--verbose', VERBOSE,
                 '--maxreads', maxreads,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE >= 2:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def make_output_folders(data_folder, adaID, VERBOSE=0):
    '''Make output folders'''
    output_filename = get_divided_filenames(data_folder, adaID, fragments=['F1'])[0]
    dirname = os.path.dirname(output_filename)
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
        if VERBOSE:
            print 'Folder created:', dirname


def test_integrity(reads, VERBOSE=True):
    '''Test integrity of read pair (should be valid at every function entry/exit)'''
    i_fwd = reads[0].is_reverse
    i_rev = not i_fwd

    # Check read integrity
    if reads[i_fwd].pos < 0:
        if VERBOSE:
            print 'Read fwd starts before 0 ('+str(reads[i_fwd].pos)+'):', reads[0].qname
        return True

    if (reads[0].mpos != reads[1].pos) or (reads[1].mpos != reads[0].pos):
        if VERBOSE:
            print 'Read pair not integer (mpos):', reads[0].qname
        return True

    if (reads[i_fwd].isize <= 0) or (reads[i_rev].isize >= 0):
        if VERBOSE:
            print 'Read pair not integer (sign of isize):', reads[0].qname
        return True
    
    if reads[i_fwd].pos + reads[i_fwd].isize != reads[i_rev].pos + \
        sum(bl for (bt, bl) in reads[i_rev].cigar if bt in (0, 2)):
        if VERBOSE:
            print 'Read pair not integer (insert size):', reads[0].qname
        return True

    if (sum(bl for (bt, bl) in reads[0].cigar if bt in (0, 1)) != reads[0].rlen) or \
       (sum(bl for (bt, bl) in reads[1].cigar if bt in (0, 1)) != reads[1].rlen):
        if VERBOSE:
            print 'Read pair not integer (CIGAR <-> seq):', reads[0].qname
        return True

    return False


def test_crossoverhang(reads, VERBOSE=True):
    '''Test reads overhanging beyond the insert size'''
    i_fwd = reads[0].is_reverse
    i_rev = not i_fwd

    if (sum(bl for (bt, bl) in reads[i_fwd].cigar if bt in (0, 2)) > reads[i_fwd].isize) or \
       (reads[i_fwd].pos > reads[i_rev].pos):
        if VERBOSE:
            print 'Read pair is cross-overhang:', reads[0].qname
        return True
    else:
        return False


def test_fragment_assignment(reads, n_frag, len_frag, VERBOSE=True):
    '''Test assignment of the read pair to a fragment'''
    # Let's assume the read is integet and has no cross-overhangs
    i_fwd = reads[0].is_reverse
    i_rev = not i_fwd

    if (reads[i_fwd].pos <  0) or \
       (reads[i_fwd].pos + reads[i_fwd].isize > len_frag):
        if VERBOSE:
            print 'Read pair is misassigned to fragment:', reads[0].qname, n_frag
        return True

    return False


def test_sanity(reads, n_frag, len_frag):
    '''Perform sanity checks on supposedly good reads'''
    # Check read integrity
    if test_integrity(reads):
        return True

    # Check cross-overhangs
    if test_crossoverhang(reads):
        return True

    # Check fragment assignment
    if test_fragment_assignment(reads, n_frag, len_frag):
        return True

    return False


def test_outer_primer(reads, pr_outs, len_reference):
    '''Test for reads starting at outer primers
    
    Note: we only test the fwd read for fwd primers, and vice versa.'''
    i_fwd = reads[0].is_reverse
    i_rev = not i_fwd

    read_fwd = reads[i_fwd]
    read_rev = reads[i_rev]

    # The first and last outer primers are special: we can spot that stuff because
    # it is mapped at 0 with an insert as first CIGAR (a read never starts with
    # a deletion) or it reaches the end of the reference with a final insert.
    if ((read_fwd.pos == 0) and (read_fwd.cigar[0][0] == 1)) or \
       ((read_fwd.pos + read_fwd.isize == len_reference) and (read_rev.cigar[-1][0] == 1)):
        return True

    # Test all fragments
    for (pr_fwd, pr_rev) in pr_outs:
        # FWD
        rfwd = np.fromstring(read_fwd.seq[:len(pr_fwd)], 'S1')
        if (rfwd == pr_fwd).mean() > 0.8:
            return True

        # REV
        rrev = np.fromstring(read_rev.seq[-len(pr_rev):], 'S1')
        if (rrev == pr_rev).mean() > 0.8:
            return True

    return False


def assign_to_fragment(reads, fragpri_pos):
    '''Assign read pair to fragments'''
    i_fwd = reads[0].is_reverse
    i_rev = not i_fwd

    # Insert coordinates
    start_fwd = reads[i_fwd].pos
    end_rev = reads[i_fwd].pos + reads[i_fwd].isize

    # If short inserts are not trimmed yet, we easily have the fwd read going
    # beyond reads[i_fwd].pos + reads[i_fwd].isize. Ignore that stuff, it's being
    # trimmed later on. For assignment, check only the insert AS A WHOLE.
    frags_pair = ((start_fwd >= fragpri_pos[0]) &
                  (start_fwd < fragpri_pos[1]) &
                  (end_rev <= fragpri_pos[1])).nonzero()[0]

    return frags_pair


def trim_crossoverhangs(reads, trim=5, include_tests=False):
    '''Trim short inserts so that they do not overhang, minus a few bases'''
    if include_tests:
        if test_integrity(reads):
            print 'trim_crossoverhangs (entry):'
            import ipdb; ipdb.set_trace()

    i_fwd = reads[0].is_reverse
    i_rev = not i_fwd

    # FWD
    end_rev = reads[i_fwd].pos + reads[i_fwd].isize
    read = reads[i_fwd]
    ref_pos = read.pos
    read_pos = 0
    cigar = []
    for i, (bt, bl) in enumerate(read.cigar):
        if bt == 0:
            if ref_pos + bl >= end_rev - trim:
                cigar.append((bt, end_rev - ref_pos - trim))
                read_pos += end_rev - ref_pos - trim
                break

            cigar.append((bt, bl))
            ref_pos += bl
            read_pos += bl

        elif bt == 1:
            cigar.append((bt, bl))
            read_pos += bl

        elif bt == 2:
            if ref_pos + bl >= end_rev - trim:
                # Do not end with a deletion (stampy would not either)
                break
            cigar.append((bt, bl))
            ref_pos += bl

        else:
            raise ValueError('CIGAR type '+str(bt)+' not recognized')
    
    seq = read.seq
    qual = read.qual
    read.seq = seq[:read_pos]
    read.qual = qual[:read_pos]
    read.cigar = cigar

    # REV (go backwards, otherwise we do not get the cigar!)
    start_fwd = reads[i_fwd].pos
    read = reads[i_rev]
    ref_pos = end_rev
    read_pos = read.rlen
    cigar = []
    for i, (bt, bl) in enumerate(read.cigar[::-1]):
        if bt == 0:
            if ref_pos - bl <= start_fwd + trim:
                cigar.append((bt, ref_pos - (start_fwd + trim)))
                read_pos -= ref_pos - (start_fwd + trim)
                ref_pos = start_fwd + trim
                break

            cigar.append((bt, bl))
            ref_pos -= bl
            read_pos -= bl

        elif bt == 1:
            cigar.append((bt, bl))
            read_pos -= bl

        elif bt == 2:
            if ref_pos - bl <= start_fwd + trim:
                # Do not end with a deletion (stampy would not either)
                break
            cigar.append((bt, bl))
            ref_pos -= bl

        else:
            raise ValueError('CIGAR type '+str(bt)+' not recognized')
    cigar = cigar[::-1]

    seq = read.seq
    qual = read.qual
    read.pos = ref_pos
    read.seq = seq[read_pos:]
    read.qual = qual[read_pos:]
    read.cigar = cigar

    # Fix mate pair
    reads[i_fwd].mpos = reads[i_rev].pos
    reads[i_rev].mpos = reads[i_fwd].pos
    isize = reads[i_rev].pos + sum(bl for bt, bl in reads[i_rev].cigar
                                   if bt in (0, 2)) - reads[i_fwd].pos
    reads[i_fwd].isize = isize
    reads[i_rev].isize = -isize

    if include_tests:
        if test_integrity(reads):
            print 'trim_crossoverhangs (exit):'
            import ipdb; ipdb.set_trace()


def trim_inner_primers(reads, frag_pos, pri, include_tests=False):
    '''Trim inner primers
    
    - frags_pos are the coordinates of the fragments trimmed of inner primers
    '''
    # Note: this function is robust against fuzzy ends, i.e. it works also if the
    # insert reads into the adapters and crap like that.
    # Because we accept short inserts, the primer can also be at the end

    if include_tests:
        if test_integrity(reads):
            print 'trim_inner_primers (entry):'
            import ipdb; ipdb.set_trace()

    tampered = False
    for read in reads:
        # FWD primer
        if read.pos < frag_pos[0]:
            tampered = True
            ref_pos = read.pos
            read_pos = 0
            cigar = read.cigar[::-1]
            for i, (bt, bl) in enumerate(read.cigar):
                if bt == 0:
                    # Strictly greater, because we dump the CIGAR ending at the
                    # right position. In the corner case, ref_pos == frag_pos[0]
                    # and we add the whole block (and do not move read_pos).
                    if ref_pos + bl > frag_pos[0]:
                        cigar[-1] = (bt, ref_pos + bl - frag_pos[0])
                        read_pos += frag_pos[0] - ref_pos
                        ref_pos = frag_pos[0]
                        break
                    cigar.pop(-1)
                    read_pos += bl
                    ref_pos += bl
                elif bt == 1:
                    cigar.pop(-1)
                    read_pos += bl
                elif bt == 2:
                    # Starting with a deletion is not allowed
                    cigar.pop(-1)
                    ref_pos += bl
                    if ref_pos > frag_pos[0]:
                        break
            cigar = cigar[::-1]

            seq = read.seq
            qual = read.qual
            read.pos = ref_pos
            read.seq = seq[read_pos:]
            read.qual = qual[read_pos:]
            read.cigar = cigar

        # REV primer
        ref_pos = read.pos + sum(bl for (bt, bl) in read.cigar if bt in (0, 2))
        if ref_pos > frag_pos[1]:
            tampered = True
            read_pos = read.rlen
            cigar = read.cigar
            for i, (bt, bl) in enumerate(read.cigar[::-1]):
                if bt == 0:
                    # Strictly less, because we dump the CIGAR starting at the
                    # right position. In the corner case, ref_pos == frag_pos[1]
                    # and we add the whole block (and do not move read_pos).
                    if ref_pos - bl < frag_pos[1]:
                        cigar[-1] = (bt, frag_pos[1] - (ref_pos - bl))
                        read_pos -= ref_pos - frag_pos[1]
                        break
                    cigar.pop(-1)
                    read_pos -= bl
                    ref_pos -= bl
                elif bt == 1:
                    cigar.pop(-1)
                    read_pos -= bl
                elif bt == 2:
                    # Ending with a deletion is not allowed
                    cigar.pop(-1)
                    ref_pos -= bl
                    if ref_pos < frag_pos[1]:
                        break

            seq = read.seq
            qual = read.qual
            read.seq = seq[:read_pos]
            read.qual = qual[:read_pos]
            read.cigar = cigar

    # Fix mate pair
    if tampered:
        i_fwd = reads[0].is_reverse
        i_rev = not i_fwd
        reads[i_fwd].mpos = reads[i_rev].pos
        reads[i_rev].mpos = reads[i_fwd].pos
        isize = reads[i_rev].pos + sum(bl for bt, bl in reads[i_rev].cigar
                                       if bt in (0, 2)) - reads[i_fwd].pos
        reads[i_fwd].isize = isize
        reads[i_rev].isize = -isize

    if include_tests:
        if test_integrity(reads):
            print 'trim_inner_primers (exit):'
            import ipdb; ipdb.set_trace()


def main_block_low_quality(reads, phred_min=20, read_len_min=50, include_tests=False,
                     VERBOSE=0):
    '''Trim low-quality stretches'''

    if include_tests:
        if test_integrity(reads):
            print 'main_block_quality (entry):'
            import ipdb; ipdb.set_trace()

    tampered = False
    for read in reads:
        # Sanger phred score used (illumina 1.8+)
        phred = np.fromstring(read.qual, np.int8) - 33
        # get all positions above the cut-off
        ind = np.asarray(phred >= phred_min, int)

        # If the whole read is below the threshold, trash the pair
        if not ind.any():
            return True

        # If the whole read is safe, do not tamper
        if ind.all():
            continue

        # We have to tamper
        tampered = True

        # divide in blocks
        switch = np.diff(ind).nonzero()[0] + 1
        ind_block_start = np.insert(switch, 0, 0)
        ind_block_end = np.append(switch, len(ind))
    
        # keep only high-q blocks
        # If the first block is high-q, even blocks are; else, odd blocks are
        first_block_good = ind[0]
        ind_block_start = ind_block_start[not first_block_good::2]
        ind_block_end = ind_block_end[not first_block_good::2]

        # get largest
        blocks_len = ind_block_end - ind_block_start
        ind_largest_block = blocks_len.argmax()

        # Check how much we lost
        if VERBOSE >= 3:
            percent_lost = 100 - 100 * (read_end - read_start) / read.rlen
            print 'Q-trim lost:', percent_lost, '%'

        # Trash tiny reads
        if blocks_len[ind_largest_block] < read_len_min:
            return True

        # rewrite read such that CIGARs are fine
        # START
        read_start = ind_block_start[ind_largest_block]
        ref_start = read.pos
        if read_start == 0:
            cigar = read.cigar
        else:
            read_pos = 0
            cigar = read.cigar[::-1]
            for (bt, bl) in read.cigar:
                # A read CAN start with an insertion
                if bt in (0, 1):
                    if read_pos + bl > read_start:
                        cigar[-1] = (bt, read_pos + bl - read_start)
                        if bt == 0:
                            ref_start += read_start - read_pos
                        break
                    cigar.pop(-1)
                    if bt == 0:
                        ref_start += bl
                    read_pos += bl
                elif bt == 2:
                    # A read cannot start with a deletion
                    cigar.pop(-1)
            cigar = cigar[::-1]

        # END (we are operating on the trimmed read now)
        read_end = ind_block_end[ind_largest_block] - read_start
        # If we go all the way, no need for trimming the end
        if read_end + read_start < read.rlen:
            read_pos = 0
            # We walk along the read this time, because it's probably faster
            # than actually trimming from the back
            cigar_new = []
            for (bt, bl) in cigar:
                # A read CAN end with an insertion
                if bt in (0, 1):
                    if read_pos + bl >= read_end:
                        cigar_new.append((bt, read_end - read_pos))
                        break
                    cigar_new.append((bt, bl))
                    read_pos += bl
                elif bt == 2:
                    # Note: a read cannot end with a deletion, so we do nothing
                    if read_pos + bl >= read_end:
                        break
                    cigar_new.append((bt, bl))

            cigar = cigar_new

        # Write properties
        seq = read.seq
        qual = read.qual
        read.pos = ref_start
        read.seq = seq[read_start: read_start + read_end]
        read.qual = qual[read_start: read_start + read_end]
        read.cigar = cigar

    # Fix mate pair stuff
    if tampered:
        i_fwd = reads[0].is_reverse
        i_rev = not i_fwd
        reads[i_fwd].mpos = reads[i_rev].pos
        reads[i_rev].mpos = reads[i_fwd].pos
        isize = reads[i_rev].pos + sum(bl for bt, bl in reads[i_rev].cigar
                                       if bt in (0, 2)) - reads[i_fwd].pos
        reads[i_fwd].isize = isize
        reads[i_rev].isize = -isize

    if include_tests:
        if test_integrity(reads):
            print 'main_block_low_quality (exit):'
            import ipdb; ipdb.set_trace()

    return False


def trim_low_quality(reads, phred_min=20, read_len_min=50, include_tests=False,
                     VERBOSE=0):
    '''Strip loq-q from left and right, but leave in the middle'''
    # The rationale of this approach is that we still have the Qs later for
    # more detailed exclusions (e.g. for minor allele frequencies)

    if include_tests:
        if test_integrity(reads):
            print 'trim_low_quality (entry):'
            import ipdb; ipdb.set_trace()

    tampered = False
    for read in reads:
        # Sanger phred score used (illumina 1.8+)
        phred = np.fromstring(read.qual, np.int8) - 33
        # get all positions above the cut-off
        ind = np.asarray(phred >= phred_min, int)

        # If the whole read is safe, do not tamper
        if ind.all():
            continue

        # Use sliding windows to check for mushy edges (we allow for ONE low-q,
        # which should cover most cases).
        read_size_min = 50
        win_size = 10
        win_qual_threshold = 9
        shift = 5
        # LEFT EDGE
        read_start = 0
        win_qual = 0
        while win_qual < win_qual_threshold:
            # If no window ever reaches the quality threshold, trash the pair
            if read_start > read.rlen - read_size_min:
                return True
            win_phred = phred[read_start: read_start + win_size]
            win_qual = (win_phred >= phred_min).sum()
            read_start += shift
        read_start -= shift
        # RIGHT EDGE
        read_end = read.rlen
        win_qual = 0
        while win_qual < win_qual_threshold:
            # If the high-q read chunk is tiny, trash the pair
            if read_end < read_start + read_size_min:
                return True
            win_phred = phred[read_end - win_size: read_end]
            win_qual = (win_phred >= phred_min).sum()
            read_end -= shift
        read_end += shift

        # If the trimmed read still has widespread low-q, it was not a trimming
        # problem: trash the pair (this happend almost never)
        if phred[read_start: read_end].mean() < 0.9:
            return True

        # If we trim nothing, proceed: this happens if the only low-q bases are
        # singletons in the middle of the read (this happens a lot, say someone
        # opened the door of the MiSeq room)
        if (read_start == 0) and (read_end == read.rlen):
            continue

        # Check how much we lost
        if VERBOSE >= 3:
            percent_lost = 100 - 100 * (read_end - read_start) / read.rlen
            print 'Q-trim lost:', percent_lost, '%'

        # or else, we have to tamper
        tampered = True

        ##FIXME
        #if read.qname == 'HWI-M01346:28:000000000-A53RP:1:1101:2303:16467':
        #    print read.qname
        #    import ipdb; ipdb.set_trace()

        # rewrite read such that CIGARs are fine
        # START
        ref_start = read.pos
        if read_start == 0:
            cigar = read.cigar
        else:
            read_pos = 0
            cigar = read.cigar[::-1]
            for (bt, bl) in read.cigar:
                # A read CAN start with an insertion
                if bt in (0, 1):
                    if read_pos + bl > read_start:
                        cigar[-1] = (bt, read_pos + bl - read_start)
                        if bt == 0:
                            ref_start += read_start - read_pos
                        break
                    cigar.pop(-1)
                    if bt == 0:
                        ref_start += bl
                    read_pos += bl
                elif bt == 2:
                    # A read cannot start with a deletion
                    cigar.pop(-1)
            cigar = cigar[::-1]

        # END (we are operating on the trimmed read now)
        read_end -= read_start
        # If we go all the way, no need for trimming the end
        if read_end + read_start < read.rlen:
            read_pos = 0
            # We walk along the read this time, because it's probably faster
            # than actually trimming from the back
            cigar_new = []
            for (bt, bl) in cigar:
                # A read CAN end with an insertion
                if bt in (0, 1):
                    if read_pos + bl >= read_end:
                        cigar_new.append((bt, read_end - read_pos))
                        break
                    cigar_new.append((bt, bl))
                    read_pos += bl
                elif bt == 2:
                    # we cannot reach read_end via a deletion, because read_pos
                    # does not increase, so there's going to be a new cigar after
                    # this (in any case, the read did never and will never end
                    # with a deletion
                    cigar_new.append((bt, bl))

            cigar = cigar_new

        # Write properties
        seq = read.seq
        qual = read.qual
        read.pos = ref_start
        read.seq = seq[read_start: read_start + read_end]
        read.qual = qual[read_start: read_start + read_end]
        read.cigar = cigar

    # Fix mate pair stuff
    if tampered:
        i_fwd = reads[0].is_reverse
        i_rev = not i_fwd
        reads[i_fwd].mpos = reads[i_rev].pos
        reads[i_rev].mpos = reads[i_fwd].pos
        isize = reads[i_rev].pos + sum(bl for bt, bl in reads[i_rev].cigar
                                       if bt in (0, 2)) - reads[i_fwd].pos

        # If extremely rare cases, we trim so much that the read becomes fully
        # cross-overhanging
        #                ------->
        #    <-----
        # we should dump the pair in this case (short inserts are dumped later anyway)
        if isize <= 0:
            return True

        reads[i_fwd].isize = isize
        reads[i_rev].isize = -isize

    if include_tests:
        if test_integrity(reads):
            print 'trim_low_quality (exit):'
            import ipdb; ipdb.set_trace()

    return False


def trim_and_divide_reads(data_folder, adaID, n_cycles, F5_primer, VERBOSE=0,
                          maxreads=-1, include_tests=False):
    '''Trim reads and divide them into fragments'''
    if VERBOSE:
        print 'Trim and divide into fragments: adaID '+'{:02d}'.format(adaID)

    # Extract fragments and primers
    fragments = ['F1', 'F2', 'F3', 'F4', F5_primer, 'F6']
    # This structure contains the fragment coordinates as of
    # - outer primers (row 0)
    # - inner primers (row 1)
    # - trimmed of all primers (row 2)
    frags_pos = np.zeros((3, 2, len(fragments)), int)
    for i, fragment in enumerate(fragments):
        pci = pcis[fragment]
        pco = pcos[fragment]
        frags_pos[0, :, i] = (pco[0][0], pco[1][1])
        frags_pos[1, :, i] = (pci[0][0], pci[1][1])
        frags_pos[2, :, i] = (pci[0][1], pci[1][0])
    # Since the reference is cropped, subtract from the positions F1 start
    # Note: now we are in the reference of the CROPPED HXB2, and start from 0!
    frags_pos -= frags_pos[1].min()
    # Make primers with masks for ambiguous nucleotides
    pr_outs = []
    for fragment in fragments:
        ptmps = primers_outer[fragment]
        for i, ptmp in enumerate(ptmps):
            ptmp = np.ma.array(list(ptmp), mask=[p not in alpha[:4] for p in ptmp])
            ptmps[i] = ptmp
        pr_outs.append(ptmps)

    # Input and output files
    input_filename = get_premapped_file(data_folder, adaID, type='bam')
    if not os.path.isfile(input_filename):
        convert_sam_to_bam(input_filename)
    output_filenames = get_divided_filenames(data_folder, adaID, fragments, type='bam')
    with pysam.Samfile(input_filename, 'rb') as bamfile:
        with pysam.Samfile(output_filenames[0], 'wb', template=bamfile) as fo_F1,\
             pysam.Samfile(output_filenames[1], 'wb', template=bamfile) as fo_F2,\
             pysam.Samfile(output_filenames[2], 'wb', template=bamfile) as fo_F3,\
             pysam.Samfile(output_filenames[3], 'wb', template=bamfile) as fo_F4,\
             pysam.Samfile(output_filenames[4], 'wb', template=bamfile) as fo_F5,\
             pysam.Samfile(output_filenames[5], 'wb', template=bamfile) as fo_F6,\
             pysam.Samfile(output_filenames[6], 'wb', template=bamfile) as fo_am,\
             pysam.Samfile(output_filenames[7], 'wb', template=bamfile) as fo_um,\
             pysam.Samfile(output_filenames[8], 'wb', template=bamfile) as fo_lq:

            # Collect the file handles
            file_handles = (fo_F1, fo_F2, fo_F3, fo_F4, fo_F5, fo_F6)

            # Iterate over the mapped reads and assign fragments
            n_mapped = [0 for fragment in fragments]
            n_unmapped = 0
            n_ambiguous = 0
            n_outer = 0
            n_lowq = 0
            for irp, reads in enumerate(pair_generator(bamfile)):

                # Stop at the maximal number of reads (for testing)
                if irp == maxreads:
                    if VERBOSE:
                        print 'Maximal number of read pairs reached:', maxreads
                    break

                # If unmapped or unpaired, mini, or insert size mini, discard
                if reads[0].is_unmapped or (not reads[0].is_proper_pair) or \
                   reads[1].is_unmapped or (not reads[1].is_proper_pair) or \
                   (reads[0].rlen < 50) or (reads[1].rlen < 50) or \
                   (np.abs(reads[0].isize) < 100):
                    if VERBOSE >= 3:
                        print 'Read pair unmapped/unpaired/tiny:', reads[0].qname
                    n_unmapped += 1
                    fo_um.write(reads[0])
                    fo_um.write(reads[1])
                    continue

                # If the insert is a misamplification from the outer primers,
                # trash it (it will have skewed amplification anyway). We cannot
                # find all of those, rather only the ones still carrying the
                # primer itself (some others have lost it while shearing). For
                # those, no matter what happens at the end (reading into adapters,
                # etc.), ONE of the reads in the pair will start exactly with one
                # outer primer: if the rev read with a rev primer, if the fwd
                # with a fwd one. Test all six.
                if test_outer_primer(reads, pr_outs, frags_pos[1, 1, -1]):
                    if VERBOSE >= 3:
                        print 'Read pair from outer primer:', reads[0].qname
                    n_outer += 1
                    fo_um.write(reads[0])
                    fo_um.write(reads[1])
                    continue

                # Assign to a fragment now, so that primer trimming is faster 
                frags_pair = assign_to_fragment(reads, frags_pos[1])

                # 1. If no fragments are possible (e.g. one read crosses the
                # fragment boundary, they map to different fragments), dump it
                if len(frags_pair) == 0:
                    n_unmapped += 1
                    fo_um.write(reads[0])
                    fo_um.write(reads[1])
                    continue

                # 2. If 2+ fragments are possible (tie), put into a special bucket
                # (essentially excluded, because we want two independent measurements
                # in the overlapping region, but we might want to recover them)
                elif len(frags_pair) > 1:
                    n_ambiguous += 1
                    fo_am.write(reads[0])
                    fo_am.write(reads[1])
                    continue

                # 3. If the intersection is a single fragment, good
                n_frag = frags_pair[0]
                trim_inner_primers(reads, frags_pos[2, :, n_frag],
                                   primers_inner[fragments[n_frag]],
                                   include_tests=include_tests)
                # Nothing bad can happen here: if the insert is small, it will have
                # at most one primer, and trimming 30 bases is fine (it is 100
                # at least!)

                # Quality trimming: if no decently long pair survives, trash
                #trashed_quality = main_block_low_quality(reads, phred_min=20, include_tests=include_tests)
                trashed_quality = trim_low_quality(reads, phred_min=20, include_tests=include_tests)
                i_fwd = reads[0].is_reverse
                if trashed_quality or (reads[i_fwd].isize < 100):
                    n_lowq += 1
                    if VERBOSE >= 3:
                        print 'Read pair has low phred quality:', reads[0].qname
                    fo_lq.write(reads[0])
                    fo_lq.write(reads[1])
                    continue

                # Check for cross-overhangs (reading into the adapters)
                #        --------------->
                #    <-----------
                # In that case, trim to perfect overlap.
                if test_crossoverhang(reads, VERBOSE=False):
                    trim_crossoverhangs(reads, trim=0, include_tests=include_tests)

                # 1. If no fragments are possible (e.g. one read crosses the
                # Change coordinates into the fragmented reference (primer-trimmed)
                for read in reads:
                    read.pos -= frags_pos[2, 0, n_frag]
                    read.mpos -= frags_pos[2, 0, n_frag]

                # Here the tests
                if include_tests:
                    if test_sanity(reads, n_frag, frags_pos[2, 1, n_frag] - frags_pos[2, 0, n_frag]):
                        print 'Tests failed:', reads[0].qname
                        import ipdb; ipdb.set_trace()

                # There we go!
                n_mapped[n_frag] += 1
                file_handles[n_frag].write(reads[0])
                file_handles[n_frag].write(reads[1])

    if VERBOSE:
        print 'Trim and divide results: adaID '+'{:02d}'.format(adaID)
        print 'Total:\t\t', irp
        print 'Mapped:\t\t', sum(n_mapped), n_mapped
        print 'Unmapped:\t', n_unmapped
        print 'Outer primer\t', n_outer
        print 'Ambiguous:\t', n_ambiguous
        print 'Low-quality:\t', n_lowq

    # Write summary to file
    with open(get_divide_summary_filename(data_folder, adaID), 'w') as f:
        f.write('Trim and divide results: adaID '+'{:02d}'.format(adaID)+'\n')
        f.write('Total:\t\t'+str(irp)+'\n')
        f.write('Mapped:\t\t'+str(sum(n_mapped))+' '+str(n_mapped)+'\n')
        f.write('Unmapped:\t'+str(n_unmapped)+'\n')
        f.write('Outer primer\t'+str(n_outer)+'\n')
        f.write('Ambiguous:\t'+str(n_ambiguous)+'\n')
        f.write('Low-quality:\t'+str(n_lowq)+'\n')


def report_coverage(data_folder, adaID, VERBOSE=0):
    '''Produce a report on rough coverage on HXB2'''
    #TODO
    pass





# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Trim and divide reads into fragments')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Maximal number of reads to analyze')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--test', action='store_true',
                        help='Include sanity checks on mapped reads (slow)')
    parser.add_argument('--report', action='store_true',
                        help='Perform quality checks and save into a report')

    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    VERBOSE = args.verbose
    maxreads = args.maxreads
    submit = args.submit
    include_tests = args.test
    report = args.report

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # Set the number of cycles of the kit (for trimming adapters in short inserts)
    n_cycles = dataset['n_cycles']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = MiSeq_runs[miseq_run]['adapters']

    # Iterate over all adaIDs
    for adaID in adaIDs:

        # Submit to the cluster self if requested
        if submit:
            if include_tests:
                raise ValueError('Tests require an interactive shell')
            fork_self(miseq_run, adaID, VERBOSE=VERBOSE, maxreads=maxreads)
            continue

        # Make output folders
        make_output_folders(data_folder, adaID, VERBOSE=VERBOSE)

        # Set the primers for fragment 5
        F5_primer = dataset['primerF5'][dataset['adapters'].index(adaID)]

        # Trim reads and assign them to a fragment (or discard)
        # Note: we pass over the file only once
        trim_and_divide_reads(data_folder, adaID, n_cycles, F5_primer,
                              maxreads=maxreads, VERBOSE=VERBOSE,
                              include_tests=include_tests)

        # Check quality and report if requested
        if report:

            # Report rough coverage on HXB2
            report_coverage(data_folder, adaID, VERBOSE=VERBOSE)

