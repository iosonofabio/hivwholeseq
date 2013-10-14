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
    from mapping.filenames import get_last_mapped
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
    from Bio.Align.Applications import MuscleCommandline
    muscle_cline = MuscleCommandline(diags=True)
    child = sp.Popen(str(muscle_cline),
                     stdin=sp.PIPE,
                     stdout=sp.PIPE,
                     stderr=sp.PIPE,
                     shell=True)
    SeqIO.write(seqs, child.stdin, "fasta")
    child.stdin.close()
    child.stderr.close()
    align = AlignIO.read(child.stdout, "fasta")
    child.stdout.close()
    return align

