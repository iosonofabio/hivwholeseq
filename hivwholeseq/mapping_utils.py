# vim: fdm=indent
'''
author:     Fabio Zanini
date:       22/08/13
content:    Settings of stampy used by our mapping scripts.
'''
# Globals
stampy_bin = '/ebio/ag-neher/share/programs/bundles/stampy-1.0.22/stampy.py'
subsrate = '0.05'
bwa_bin = '/ebio/ag-neher/share/programs/bin/bwa'
spades_bin = '/ebio/ag-neher/share/programs/bundles/SPAdes-2.5.0-Linux/bin/spades.py'




# Functions
def pair_generator(iterable):
    '''Generator for pairs in interleaved files'''
    # Note: the last item is lost if odd
    it = iter(iterable)
    while True:
        try:
            p1 = it.next()
            p2 = it.next()
            yield (p1, p2)
        except StopIteration:
            raise


def get_ind_good_cigars(cigar, match_len_min=20, full_output=False):
    '''Keep only CIGAR blocks between two long matches'''
    from numpy import array

    criterion = lambda x: (x[0] == 0) and (x[1] >= match_len_min)
    good_cigars = array(map(criterion, cigar), bool, ndmin=1)

    # If there are 2+ good CIGARs, keep also stuff in between
    n_good_cigars = (good_cigars).sum()
    if n_good_cigars >= 2:
        tmp = good_cigars.nonzero()[0]
        first_good_cigar = tmp[0]
        last_good_cigar = tmp[-1]
        good_cigars[first_good_cigar: last_good_cigar + 1] = True
    elif n_good_cigars == 1:
        first_good_cigar = last_good_cigar = good_cigars.nonzero()[0][0]
    else:
        first_good_cigar = last_good_cigar = None

    if full_output:
        return good_cigars, first_good_cigar, last_good_cigar
    else:
        return good_cigars


def get_range_good_cigars(cigar, pos, match_len_min=30,
                          trim_left=3, trim_right=3):
    '''Get the range of the good CIGARs, in the read and ref sys of coordinates'''
    from numpy import array

    criterion = lambda x: (x[0] == 0) and (x[1] >= match_len_min)
    good_cigars = array(map(criterion, cigar), bool, ndmin=1)
    if not (good_cigars).any():
        return (None, None)
    else:
        tmp = good_cigars.nonzero()[0]
        first_good_cigar = tmp[0]
        last_good_cigar = tmp[-1]

        # Get the start
        start_read = 0
        start_ref = pos	# pysam already switches to the 0-start system
        for (block_type, block_len) in cigar[:first_good_cigar]:
            if block_type == 0:
                start_read += block_len
                start_ref += block_len
            elif block_type == 1:
                start_read += block_len
            elif block_type == 2:
                start_ref += block_len
            else:
                raise ValueError('CIGAR type '+str(block_type)+' not recognized')

        # Get the end
        end_read = start_read
        end_ref = start_ref
        for (block_type, block_len) in cigar[first_good_cigar: last_good_cigar + 1]:
            if block_type == 0:
                end_read += block_len
                end_ref += block_len
            elif block_type == 1:
                end_read += block_len
            elif block_type == 2:
                end_ref += block_len
            else:
                raise ValueError('CIGAR type '+str(block_type)+' not recognized')

        # If some CIGARs are chopped off, trim a few bases from that side too
        # This is also useful to avoid short random matches to HXB2 that appear
        # like crossing the fragment boundary but are actually only a short
        # insert reading back into the illumina adapters
        # Note also that trimming is fine in both coordinate systems because at
        # the edges there is a long match block:
        #    match_len_min > trim_left + trim_right
        # must be fulfilled, otherwise we do bogus.
        if first_good_cigar != 0:
            start_read += trim_left
            start_ref += trim_left
        if last_good_cigar != len(cigar) - 1:
            end_read -= trim_right
            end_ref -= trim_right

        return ((start_read, end_read), (start_ref, end_ref))


def get_trims_from_good_cigars(good_cigars, trim_left=3, trim_right=3):
    '''Set the trimming of cigars'''
    from numpy import zeros

    trims_left_right = zeros((len(good_cigars), 2), int)

    if good_cigars.any():
        tmp = good_cigars.nonzero()[0] 
        first_good_cigar = tmp[0]
        last_good_cigar = tmp[-1]
        if first_good_cigar != 0:
            trims_left_right[first_good_cigar, 0] = trim_left
        if last_good_cigar != len(good_cigars) - 1:
            trims_left_right[last_good_cigar, 1] = trim_right

    return trims_left_right


def convert_sam_to_bam(bamfilename, samfilename=None):
    '''Convert SAM file to BAM file format'''
    import pysam
    if samfilename is None:
        samfilename = bamfilename[:-3]+'sam'

    samfile = pysam.Samfile(samfilename, 'r')
    bamfile = pysam.Samfile(bamfilename, 'wb', template=samfile)
    for s in samfile: bamfile.write(s)
    samfile.close()
    bamfile.close()


def convert_bam_to_sam(samfilename, bamfilename=None):
    '''Convert BAM file to SAM file format'''
    import pysam
    if bamfilename is None:
        bamfilename = samfilename[:-3]+'bam'

    bamfile = pysam.Samfile(bamfilename, 'rb')
    samfile = pysam.Samfile(samfilename, 'w', template=bamfile)
    for s in bamfile: samfile.write(s)
    bamfile.close()
    samfile.close()


def get_fragment_list(data_folder, adaID):
    '''Get the sorted list of fragments as of the BAM file'''
    import pysam
    from hivwholegenome.filenames import get_last_mapped
    bamfilename = get_last_mapped(data_folder, adaID)
    bamfile = pysam.Samfile(bamfilename, 'rb')
    chromosomes = bamfile.references
    bamfile.close()
    return chromosomes


def get_read_start_end(read):
    '''Get the start and end position of a read in its reference'''
    start = read.pos
    len_ref = sum(block_len for (block_type, block_len) in read.cigar
                  if block_type in [0, 2])
    end = start + len_ref
    return (start, end)


def reads_to_seqrecord(reads):
    '''Build a FASTQ record out of BAM reads
    
    Note: copied from Bio.SeqIO.QualityIO.py
    '''
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    # Precompute conversion table
    SANGER_SCORE_OFFSET = ord("!")
    q_mapping = dict()
    for letter in range(0, 255):
        q_mapping[chr(letter)] = letter - SANGER_SCORE_OFFSET
    
    seqs = []
    for read in reads:
        # Get the sequence first
        descr = read.qname
        id = read.qname
        name = id
        from Bio.Alphabet import IUPAC
        record = SeqRecord(Seq(read.seq, IUPAC.ambiguous_dna),
                           id=id, name=name, description=descr)
    
        # Get the qualities second
        qualities = [q_mapping[letter] for letter in read.qual]
        if qualities and (min(qualities) < 0 or max(qualities) > 93):
            raise ValueError("Invalid character in quality string")

        #For speed, will now use a dirty trick to speed up assigning the
        #qualities. We do this to bypass the length check imposed by the
        #per-letter-annotations restricted dict (as this has already been
        #checked by FastqGeneralIterator). This is equivalent to:
        #record.letter_annotations["phred_quality"] = qualities
        dict.__setitem__(record._per_letter_annotations,
                         "phred_quality", qualities)

        seqs.append(record)

    return seqs


def sort_bam(bamfilename_sorted, bamfilename_unsorted=None):
    '''Sort BAM file'''
    import pysam

    if bamfilename_unsorted is None:
        bamfilename_unsorted = bamfilename_sorted[:-11]+'.bam'

    pysam.sort(bamfilename_unsorted, bamfilename_sorted[:-4])


def index_bam(bamfilename_sorted):
    '''Index a BAM file'''
    import pysam

    pysam.index(bamfilename_sorted)


def align_muscle(*seqs):
    '''Global alignment of sequences via MUSCLE'''
    import subprocess as sp
    from Bio import AlignIO, SeqIO
    from Bio.Align.Applications import MuscleCommandline
    muscle_cline = MuscleCommandline(diags=True, quiet=True)
    child = sp.Popen(str(muscle_cline),
                     stdin=sp.PIPE,
                     stdout=sp.PIPE,
                     stderr=sp.PIPE,
                     shell=True)
    SeqIO.write(seqs, child.stdin, "fasta")
    child.stdin.close()
    align = AlignIO.read(child.stdout, "fasta")
    child.stderr.close()
    child.stdout.close()
    return align


def get_number_reads(bamfilename, format='bam'):
    '''Count the reads (not pairs) in a BAM/SAM file'''
    import pysam
    file_modes = {'bam': 'rb', 'sam': 'r'}
    with pysam.Samfile(bamfilename, file_modes[format]) as bamfile:
        n_reads = sum(1 for read in bamfile)
    return n_reads


def test_read_pair_integrity(read_pair, VERBOSE=True):
    '''Test integrity of read pair'''
    i_fwd = read_pair[0].is_reverse
    i_rev = not i_fwd

    if read_pair[i_fwd].pos < 0:
        if VERBOSE:
            print 'Read fwd starts before 0 ('+str(read_pair[i_fwd].pos)+'):', read_pair[0].qname
        return True

    if (read_pair[0].mpos != read_pair[1].pos) or (read_pair[1].mpos != read_pair[0].pos):
        if VERBOSE:
            print 'Read pair not integer (mpos):', read_pair[0].qname
        return True

    if (read_pair[i_fwd].isize <= 0) or (read_pair[i_rev].isize >= 0):
        if VERBOSE:
            print 'Read pair not integer (sign of isize):', read_pair[0].qname
        return True
    
    if read_pair[i_fwd].pos + read_pair[i_fwd].isize != read_pair[i_rev].pos + \
        sum(bl for (bt, bl) in read_pair[i_rev].cigar if bt in (0, 2)):
        if VERBOSE:
            print 'Read pair not integer (insert size):', read_pair[0].qname
        return True

    if (sum(bl for (bt, bl) in read_pair[0].cigar if bt in (0, 1)) != read_pair[0].rlen) or \
       (sum(bl for (bt, bl) in read_pair[1].cigar if bt in (0, 1)) != read_pair[1].rlen):
        if VERBOSE:
            print 'Read pair not integer (CIGAR <-> seq):', read_pair[0].qname
        return True

    return False


def test_read_pair_crossoverhang(read_pair, VERBOSE=True):
    '''Test read_pair overhanging beyond the insert size'''
    i_fwd = read_pair[0].is_reverse
    i_rev = not i_fwd

    if (sum(bl for (bt, bl) in read_pair[i_fwd].cigar if bt in (0, 2)) > read_pair[i_fwd].isize) or \
       (read_pair[i_fwd].pos > read_pair[i_rev].pos):
        if VERBOSE:
            print 'Read pair is cross-overhang:', read_pair[0].qname
        return True
    else:
        return False


def main_block_read_pair_low_quality(read_pair, phred_min=20, read_len_min=50, include_tests=False,
                     VERBOSE=0):
    '''Keep only the largest high-phred block of each read in an insert'''
    import numpy as np

    if include_tests:
        if test_read_pair_integrity(read_pair):
            print 'main_block_quality (entry):'
            import ipdb; ipdb.set_trace()

    tampered = False
    for read in read_pair:
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

        # Trash tiny read_pair
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
        i_fwd = read_pair[0].is_reverse
        i_rev = not i_fwd
        read_pair[i_fwd].mpos = read_pair[i_rev].pos
        read_pair[i_rev].mpos = read_pair[i_fwd].pos
        isize = read_pair[i_rev].pos + \
                sum(bl for bt, bl in read_pair[i_rev].cigar if bt in (0, 2)) - \
                read_pair[i_fwd].pos
        read_pair[i_fwd].isize = isize
        read_pair[i_rev].isize = -isize

    if include_tests:
        if test_read_pair_integrity(read_pair):
            print 'main_block_low_quality (exit):'
            import ipdb; ipdb.set_trace()

    return False


def trim_read_pair_low_quality(read_pair,
                               phred_min=20,
                               read_len_min=50,
                               include_tests=False,
                               VERBOSE=0):
    '''Strip loq-phred from left and right edges, but leave in the middle'''
    # The rationale of this approach is that we still have the Qs later for
    # more detailed exclusions (e.g. for minor allele frequencies)
    import numpy as np

    if include_tests:
        if test_read_pair_integrity(read_pair):
            print 'trim_low_quality (entry):'
            import ipdb; ipdb.set_trace()

    tampered = False
    for read in read_pair:
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
        # NOTE: changed recently (20/01/14)
        if (phred[read_start: read_end] >= phred_min).mean() < 0.9:
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
        i_fwd = read_pair[0].is_reverse
        i_rev = not i_fwd
        read_pair[i_fwd].mpos = read_pair[i_rev].pos
        read_pair[i_rev].mpos = read_pair[i_fwd].pos
        isize = read_pair[i_rev].pos + \
                sum(bl for bt, bl in read_pair[i_rev].cigar if bt in (0, 2)) - \
                read_pair[i_fwd].pos

        # If extremely rare cases, we trim so much that the read becomes fully
        # cross-overhanging
        #                ------->
        #    <-----
        # we should dump the pair in this case (short inserts are dumped later anyway)
        if isize <= 0:
            return True

        read_pair[i_fwd].isize = isize
        read_pair[i_rev].isize = -isize

    if include_tests:
        if test_read_pair_integrity(read_pair):
            print 'trim_low_quality (exit):'
            import ipdb; ipdb.set_trace()

    return False


def trim_read_pair_crossoverhangs(read_pair, trim=5, include_tests=False):
    '''Trim short inserts so that they do not overhang, minus a few bases'''
    if include_tests:
        if test_read_pair_integrity(read_pair):
            print 'trim_crossoverhangs (entry):'
            import ipdb; ipdb.set_trace()

    i_fwd = read_pair[0].is_reverse
    i_rev = not i_fwd

    # FWD
    end_rev = read_pair[i_fwd].pos + read_pair[i_fwd].isize
    read = read_pair[i_fwd]
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
    start_fwd = read_pair[i_fwd].pos
    read = read_pair[i_rev]
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
    read_pair[i_fwd].mpos = read_pair[i_rev].pos
    read_pair[i_rev].mpos = read_pair[i_fwd].pos
    isize = read_pair[i_rev].pos + \
            sum(bl for bt, bl in read_pair[i_rev].cigar if bt in (0, 2)) -\
            read_pair[i_fwd].pos
    read_pair[i_fwd].isize = isize
    read_pair[i_rev].isize = -isize

    if include_tests:
        if test_read_pair_integrity(read_pair):
            print 'trim_crossoverhangs (exit):'
            import ipdb; ipdb.set_trace()



