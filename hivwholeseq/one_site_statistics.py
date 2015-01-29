# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/11/13
content:    Collection of functions to do single site statistics (allele counts,
            coverage, allele frequencies).
'''
# Modules
from collections import defaultdict, Counter
import numpy as np
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import ambiguous_dna

from .utils.sequence import alpha, alphaa
from .miseq import read_types
from .utils.mapping import get_ind_good_cigars
from .utils.mapping import align_muscle


# Functions
def get_allele_counts_read(read, counts_out, inserts_out,
                           qual_min=30, length=None, VERBOSE=0):
    '''Get allele counts and insertions from a single read

    Parameters:
       counts_out (ndarray, alphabet x sequence length): output data structure for counts
       inserts_out (nested defaultdict): output data structure for insertions
    '''

    # Read CIGARs (they should be clean by now)
    seq = np.fromstring(read.seq, 'S1')
    qual = np.fromstring(read.qual, np.int8) - 33
    pos = read.pos

    # Iterate over CIGARs
    for ic, (block_type, block_len) in enumerate(read.cigar):

        # Check for pos: it should never exceed the length of the fragment
        if (length is not None) and (block_type in [0, 1, 2]) and (pos >= length):
            raise ValueError('Pos exceeded the length of the fragment')
    
        # Inline block
        if block_type == 0:
            seqb = seq[:block_len]
            qualb = qual[:block_len]
            # Increment counts
            for j, a in enumerate(alpha):
                posa = ((seqb == a) & (qualb >= qual_min)).nonzero()[0]
                if len(posa):
                    counts_out[j, pos + posa] += 1
    
            # Chop off this block
            if ic != len(read.cigar) - 1:
                seq = seq[block_len:]
                qual = qual[block_len:]
                pos += block_len
    
        # Deletion
        elif block_type == 2:
            # Increment gap counts
            counts_out[4, pos:pos + block_len] += 1
    
            # Chop off pos, but not sequence
            pos += block_len
    
        # Insertion
        # an insert @ pos 391 means that seq[:391] is BEFORE the insert,
        # THEN the insert, FINALLY comes seq[391:]
        elif block_type == 1:
            seqb = seq[:block_len]
            qualb = qual[:block_len]
            # Accept only high-quality inserts
            if (qualb >= qual_min).all():
                inserts_out[pos][seqb.tostring()] += 1
    
            # Chop off seq, but not pos
            if ic != len(read.cigar) - 1:
                seq = seq[block_len:]
                qual = qual[block_len:]
    
        # Other types of cigar?
        else:
            raise ValueError('CIGAR type '+str(block_type)+' not recognized')


def get_allele_counts_aa_read(read, start, end, counts_out, qual_min=30,
                              VERBOSE=0):
    '''Get allele counts as amino acids from a single read.

    NOTE: It is assumed that (end - start) has a length of a multiple of 3.
    This function does not check for it again.

    Parameters:
       counts_out (ndarray, alphabet x protein length): output data structure for counts
    '''
    from Bio.Seq import translate

    # Read CIGARs (they should be clean by now)
    seq = np.fromstring(read.seq, 'S1')
    qual = np.fromstring(read.qual, np.int8) - 33
    pos = read.pos
    pos_end = pos + sum(bl for (bt, bl) in read.cigar if bt in (0, 2))

    # If the read does not cover, skip
    if (pos > end - 2) or (pos_end < start):
        return

    # Iterate over CIGARs
    pos_ref = pos
    pos_read = 0
    for ic, (bt, bl) in enumerate(read.cigar):
        if bt == 1:
            pos_read += bl
            continue

        # 1. we have not reached the start, just move on
        if pos_ref + bl < start:
            pos_ref += bl
            if bt == 0:
                pos_read += bl

            # Check we are not beyond the end
            if pos_ref > end - 2:
                return

            continue

        # 2. we enter the protein now
        if pos_ref <= start:
            startb = start - pos_ref

        # 3. we are already in the protein
        else:
            # Cut block to full codons from the left
            # NOTE: modulo works for negative numbers in Python
            startb = (start - pos_ref) % 3
        

        # Cut block to full codons from the right
        endb = pos_ref + bl
        endb -= (endb - startb) % 3
        lb = endb - startb

        # If the block is less than one amino acid, skip
        if lb < 3:
            pos_ref += bl
            if bt == 0:
                pos_read += bl

            # Check we are not beyond the end
            if pos_ref > end - 2:
                return

        # Get output amino acid coordinates for the block
        start_pr = ((pos_ref + startb) - start) // 3
        end_pr = start_pr + ((endb - startb) // 3)

        if bt == 2:
            # We assume gaps is at the second-last position
            counts_out[-2, start_pr: end_pr] += 1

        else:
            # Ask for minimal quality at all three codon positions
            qualb = (qual[pos_read + startb: pos_read + endb]
                     .reshape((lb // 3, 3))
                     .min(axis=1))

            seqb = seq[pos_read + startb: pos_read + endb]
            aab = np.fromstring(translate(''.join(seqb)), 'S1')

            for j, a in enumerate(alphaa):
                posa = ((aab == a) & (qualb >= qual_min)).nonzero()[0]
                if len(posa):
                    counts_out[j, start_pr + posa] += 1

            pos_read += bl

        pos_ref += bl

        # Check we are not beyond the end
        if pos_ref > end - 2:
            return





def get_allele_counts_insertions_from_file(bamfilename, length, qual_min=30,
                                           maxreads=-1, VERBOSE=0):
    '''Get the allele counts and insertions'''
    from collections import defaultdict, Counter

    # Prepare output structures
    counts = np.zeros((len(read_types), len(alpha), length), int)
    # Note: the data structure for inserts is a nested dict with:
    # read type --> position --> string --> count
    #  (list)        (dict)      (dict)     (int)
    inserts = [defaultdict(lambda: Counter()) for rt in read_types]

    # Open BAM file
    # Note: the reads should already be filtered of unmapped stuff at this point
    with pysam.Samfile(bamfilename, 'rb') as bamfile:

        # Iterate over single reads
        for i, read in enumerate(bamfile):

            # Max number of reads
            if i == maxreads:
                if VERBOSE >= 2:
                    print 'Max reads reached:', maxreads
                break
        
            # Print output
            if (VERBOSE >= 3) and (not ((i +1) % 10000)):
                print (i+1)
        
            # Divide by read 1/2 and forward/reverse
            js = 2 * read.is_read2 + read.is_reverse

            get_allele_counts_read(read, counts[js], inserts[js],
                                   length=length,
                                   qual_min=qual_min,
                                   VERBOSE=VERBOSE)

    return (counts, inserts)


def get_allele_counts_insertions_from_file_unfiltered(bamfilename, length, qual_min=30,
                                                      match_len_min=10,
                                                      maxreads=-1, VERBOSE=0):
    '''Get the allele counts and insertions
    
    Parameters:
       - maxreads: limit the counts to a random subset of the reads of this size
    '''
    # Prepare output structures
    counts = np.zeros((len(read_types), len(alpha), length), int)
    # Note: the data structure for inserts is a nested dict with:
    # position --> string --> read type --> count
    #  (dict)      (dict)       (list)      (int)
    inserts = defaultdict(lambda: defaultdict(lambda: np.zeros(len(read_types), int)))

    # Open BAM file
    # Note: the reads should already be filtered of unmapped stuff at this point
    with pysam.Samfile(bamfilename, 'rb') as bamfile:

        if maxreads != -1:
            from hivwholeseq.utils.mapping import extract_mapped_reads_subsample_open
            read_iter = extract_mapped_reads_subsample_open(bamfile, maxreads,
                                                            VERBOSE=VERBOSE,
                                                            pairs=False)
        else:
            read_iter = bamfile

        # Iterate over single reads
        for i, read in enumerate(read_iter):

            # Max number of reads
            if i == maxreads:
                if VERBOSE >= 2:
                    print 'Max reads reached:', maxreads
                break
        
            # Print output
            if (VERBOSE >= 3) and (not ((i +1) % 10000)):
                print (i+1)

            # NOTE: since we change the consensus all the time, mapping is never
            # safe, and we have to filter the results thoroughly.

            # If unmapped/unpaired, trash
            if read.is_unmapped or (not read.is_proper_pair) or (read.isize == 0):
                if VERBOSE >= 3:
                        print 'Read '+read.qname+': unmapped/unpaired/no isize'
                continue

            # Get good CIGARs
            (good_cigars, first_good_cigar, last_good_cigar) = \
                    get_ind_good_cigars(read.cigar, match_len_min=match_len_min,
                                        full_output=True)

            # If no good CIGARs, give up
            if not good_cigars.any():
                continue
                    
            # Divide by read 1/2 and forward/reverse
            js = 2 * read.is_read2 + read.is_reverse
        
            # Read CIGARs
            seq = np.fromstring(read.seq, 'S1')
            qual = np.fromstring(read.qual, np.int8) - 33
            pos = read.pos
            cigar = read.cigar
            len_cig = len(cigar)            

            # Iterate over CIGARs
            for ic, (block_type, block_len) in enumerate(cigar):

                # Check for pos: it should never exceed the length of the fragment
                if (block_type in [0, 1, 2]) and (pos > length):
                    raise ValueError('Pos exceeded the length of the fragment')
            
                # Inline block
                if block_type == 0:
                    # Keep only stuff from good CIGARs
                    if first_good_cigar <= ic <= last_good_cigar:
                        seqb = seq[:block_len]
                        qualb = qual[:block_len]
                        # Increment counts
                        for j, a in enumerate(alpha):
                            posa = ((seqb == a) & (qualb >= qual_min)).nonzero()[0]
                            if len(posa):
                                counts[js, j, pos + posa] += 1
            
                    # Chop off this block
                    if ic != len_cig - 1:
                        seq = seq[block_len:]
                        qual = qual[block_len:]
                        pos += block_len
            
                # Deletion
                elif block_type == 2:
                    # Keep only stuff from good CIGARs
                    if first_good_cigar <= ic <= last_good_cigar:

                        # Increment gap counts
                        counts[js, 4, pos:pos + block_len] += 1
            
                    # Chop off pos, but not sequence
                    pos += block_len
            
                # Insertion
                # an insert @ pos 391 means that seq[:391] is BEFORE the insert,
                # THEN the insert, FINALLY comes seq[391:]
                elif block_type == 1:
                    # Keep only stuff from good CIGARs
                    if first_good_cigar <= ic <= last_good_cigar:
                        seqb = seq[:block_len]
                        qualb = qual[:block_len]
                        # Accept only high-quality inserts
                        if (qualb >= qual_min).all():
                            inserts[pos][seqb.tostring()][js] += 1
            
                    # Chop off seq, but not pos
                    if ic != len_cig - 1:
                        seq = seq[block_len:]
                        qual = qual[block_len:]
            
                # Other types of cigar?
                else:
                    raise ValueError('CIGAR type '+str(block_type)+' not recognized')

    return (counts, inserts)


def filter_nus(counts, coverage=None, VERBOSE=0):
    '''Filter allele frequencies from the four read types'''
    from scipy.stats import chi2_contingency

    if coverage is None:
        coverage = counts.sum(axis=1)

    # Divide binarily
    nocounts = (coverage - counts.swapaxes(0, 1)).swapaxes(0, 1)

    # Set counts and similia: sum read1 and read2
    counts_f = counts[0] + counts[2]
    counts_b = counts[1] + counts[3]
    nocounts_f = nocounts[0] + nocounts[2]
    nocounts_b = nocounts[1] + nocounts[3]
    cov_f = coverage[0] + coverage[2]
    cov_b = coverage[1] + coverage[3]
    ind_low_cov_f = cov_f < 10
    ind_low_cov_b = cov_b < 10
    ind_high_cov_both = (-ind_low_cov_f) & (-ind_low_cov_b)

    nu_filtered = np.ma.masked_all((len(alpha), counts.shape[-1]))

    # Iterate over positions
    for i in xrange(counts.shape[-1]):
        
        # 1. if we cover neither fwd nor rev, keep masked
        if ind_low_cov_f[i] and ind_low_cov_b[i]:
            if VERBOSE >= 4:
                print 'pos', i, 'not covered'
            pass

        # 2. if we cover only one of them, well, just take the
        # arithmetic sum of counts
        elif ind_low_cov_f[i] != ind_low_cov_b[i]:
            nu_filtered[:, i] = 1.0 * counts[:, :, i].sum(axis=0) / coverage[:, i].sum()
            if VERBOSE >= 4:
                print 'pos', i, 'covered only once'
                
        # 3. If we cover both, check whether the counts are significantly different
        else:

            # Check all alleles
            for j in xrange(len(alpha)):
                # To make a table, you must have coverage for both
                cm = np.array([[counts_f[j, i], nocounts_f[j, i]],
                               [counts_b[j, i], nocounts_b[j, i]]], int)
                chi2, pval = chi2_contingency(cm + 1)[:2]
    
                # If they are not significantly different, sum the counts
                if (pval > 1e-6):
                    nu_filtered[j, i] = 1.0 * counts[:, j, i].sum(axis=0) / coverage[:, i].sum()
                # If they are different by a significant and reasonable amount, take
                # the value further away from 0.5
                else:
                    nu_f = 1.0 * counts_f[j, i] / cov_f[i]
                    nu_b = 1.0 * counts_b[j, i] / cov_b[i]
                    if np.abs(nu_f - 0.5) > np.abs(nu_b - 0.5):
                        nu_filtered[j, i] = nu_f
                    else:
                        nu_filtered[j, i] = nu_b

                    if VERBOSE >= 3:
                        print 'pos', i, 'base', alpha[j], 'nu_f', nu_f, 'nu_b', nu_b

    # Renormalize to 1
    nu_filtered /= nu_filtered.sum(axis=0)

    # Get rid of the mask if not needed
    if not nu_filtered.mask.any():
        nu_filtered = nu_filtered.data

    return nu_filtered


def get_minor_allele_frequencies(afs, alpha=None):
    '''Get the identity and frequency of the top minor allele at every site'''
    if alpha is None:
        from hivwholeseq.miseq import alpha

    allm = np.zeros(afs.shape[1], 'S1')
    num = np.zeros(afs.shape[1])

    for i, af in enumerate(afs.T):
        ind = np.argsort(af)[-2]
        allm = alpha[ind]
        num[i] = af[ind]
    return allm, num


def build_consensus_from_allele_counts_insertions(counts, inserts, 
                                                  coverage_min=10,
                                                  align_inserts=False,
                                                  VERBOSE=0):
    '''Build consensus sequence (including insertions)
    
    Parameters:
        - coverage_min: minimal coverage required to not be considered N
        - align_inserts: make a MSA of insers at each position intead of using prefix/suffix trees
    '''
    import re
    from collections import Counter
    from operator import itemgetter

    # Make allele count consensi for each of the four categories (fwd/rev, r1/r2)
    consensi = np.zeros((counts.shape[0], counts.shape[-1]), 'S1')
    for js, count in enumerate(counts):
        # Positions without reads are considered N
        # (this should happen only at the ends)
        count.T[(count).sum(axis=0) == 0] = np.array([0, 0, 0, 0, 0, 1])
        consensi[js] = alpha[count.argmax(axis=0)]

    # Make final consensus
    # This has two steps: 1. counts; 2. inserts
    # We should check that different read types agree (e.g. forward/reverse) on
    # both counts and inserts if they are covered
    # 1.0: Prepare everything as 'N': ambiguous
    consensus = np.repeat('N', counts.shape[-1])
    # 1.1: Put safe stuff
    ind_agree = (consensi == consensi[0]).all(axis=0)
    consensus[ind_agree] = consensi[0, ind_agree]
    # 1.2: Ambiguous stuff requires more care
    polymorphic = []
    amb_pos = (-ind_agree).nonzero()[0]
    for pos in amb_pos:
        cons_pos = consensi[:, pos]
        # Limit to read types that have information: if none have info, we
        # have done this already (read types agree, albeit on no info)
        cons_pos = cons_pos[cons_pos != 'N']

        # If there is unanimity, ok
        if (cons_pos == cons_pos[0]).all():
            consensus[pos] = cons_pos[0]
        # else, assign only likely mismapped deletions
        else:
            # Restrict to nongapped things (caused by bad mapping)
            cons_pos = cons_pos[cons_pos != '-']
            if (cons_pos == cons_pos[0]).all():
                consensus[pos] = cons_pos[0]
            # In case of polymorphisms, take any most abundant nucleotide
            else:
                polymorphic.append(pos)
                tmp = zip(*Counter(cons_pos).items())
                consensus[pos] = tmp[0][np.argmax(tmp[1])]
    
    # 2. Inserts
    # This is a bit tricky, because we could have a variety of inserts at the
    # same position because of PCR or sequencing errors (what we do with that
    # stuff, is another question). So, assuming the mapping program maps
    # everything decently at the same starting position, we have to look for
    # prefix families
    insert_consensus = []

    # Iterate over all insert positions
    for pos, insert in inserts.iteritems():

        # Get the coverage around the insert in all four read types
        # Note: the sum is over the four nucleotides, the mean over the two positions
        cov = counts[:, :, pos: min(counts.shape[-1], pos + 2)].sum(axis=1).mean(axis=1)

        # Determine what read types have coverage
        is_cov = cov > coverage_min
        covcs = cov[is_cov]

        # A single covered read type is sufficient: we are going to do little
        # with that, but as a consensus it's fine
        if not is_cov.any():
            continue

        if align_inserts:
            # Take the ABUNDANT insertions at this positions
            insert_big = [key for (key, val) in insert.iteritems() \
                          if (val[is_cov] > 0.1 * covcs).any()]
            if not insert_big:
                continue

            # Get abundances
            val_big = np.array([insert[key] for key in insert_big], int)

            # Make MSA
            ali = align_muscle(*[SeqRecord(Seq(key, ambiguous_dna), id=str(i)) \
                                 for (i, key) in enumerate(insert_big)], sort=True)
            ali = np.array(ali, 'S1')

            # Get allele counts
            cts_ins = np.zeros((len(read_types), len(alpha), ali.shape[-1]), int)
            for pos_ali in xrange(ali.shape[-1]):
                for j, a in enumerate(alpha):
                    ind_ins = (ali[:, pos_ali] == a).nonzero()[0]
                    if len(ind_ins):
                        cts_ins[:, j, pos_ali] = val_big[ind_ins].sum(axis=0)

            # Sum read types
            cts_ins = cts_ins.sum(axis=0)

            # Get major allele counts
            cts_ins_maj = np.max(cts_ins, axis=0)
            all_ins_maj = np.argmax(cts_ins, axis=0)

            # Get only the ones above 50%
            ins = ''.join(alpha[all_ins_maj[cts_ins_maj > 0.5 * cov.sum()]])

            if ins:
                insert_consensus.append((pos, ins))

        else:

            # Use prefixes/suffixes, and come down from highest frequencies
            # (averaging fractions is not a great idea, but -- oh, well).
            # Implement it as a count table: this is not very efficient a priori,
            # but Python is lame anyway if we start making nested loops
            len_max = max(map(len, insert.keys()))
            align_pre = np.tile('-', (len(insert), len_max))
            align_suf = np.tile('-', (len(insert), len_max))
            counts_loc = np.zeros((len(insert), len(read_types)), int)
            for i, (s, cs) in enumerate(insert.iteritems()):
                align_pre[i, :len(s)] = list(s)
                align_suf[i, -len(s):] = list(s)
                counts_loc[i] = cs
    
            # Make the allele count tables
            counts_pre_table = np.zeros((len(read_types), len(alpha), len_max), int)
            counts_suf_table = np.zeros((len(read_types), len(alpha), len_max), int)
            for j, a in enumerate(alpha):
                counts_pre_table[:, j] = np.dot(counts_loc.T, (align_pre == a))
                counts_suf_table[:, j] = np.dot(counts_loc.T, (align_suf == a))
    
            # Look whether any position of the insertion has more than 50% freq
            # The prefix goes: ----> x
            ins_pre = []
            for cs in counts_pre_table.swapaxes(0, 2):
                freq = 1.0 * (cs[:, is_cov] / covcs).mean(axis=1)
                iaM = freq.argmax()
                if (alpha[iaM] != '-') and (freq[iaM] > 0.5):
                    ins_pre.append(alpha[iaM])
                else:
                    break
            # The suffix goes: x <----
            ins_suf = []
            for cs in counts_suf_table.swapaxes(0, 2)[::-1]:
                freq = 1.0 * (cs[:, is_cov] / covcs).mean(axis=1)
                iaM = freq.argmax()
                if (alpha[iaM] != '-') and (freq[iaM] > 0.5):
                    ins_suf.append(alpha[iaM])
                else:
                    break
            ins_suf.reverse()
    
            if VERBOSE >= 4:
                if ins_pre or ins_suf:
                    print ''.join(ins_pre)
                    print ''.join(ins_suf)
    
            # Add the insertion to the list (they should agree in most cases)
            if ins_pre:
                insert_consensus.append((pos, ''.join(ins_pre)))
            elif ins_suf:
                insert_consensus.append((pos, ''.join(ins_suf)))

    # 3. put inserts in
    insert_consensus.sort(key=itemgetter(0))
    consensus_final = []
    pos = 0
    for insert_name in insert_consensus:
        # Indices should be fine...
        # an insert @ pos 391 means that seq[:391] is BEFORE the insert,
        # THEN the insert, FINALLY comes seq[391:]
        consensus_final.append(''.join(consensus[pos:insert_name[0]]))
        consensus_final.append(insert_name[1])
        pos = insert_name[0]
    consensus_final.append(''.join(consensus[pos:]))
    # Strip initial and final Ns and gaps
    consensus_final = ''.join(consensus_final).strip('N')
    consensus_final = re.sub('-', '', consensus_final)

    return consensus_final


def build_consensus_from_mapped_reads(bamfilename, maxreads=2000, block_len=100,
                                      min_reads_per_group=30,
                                      len_fragment=2000, VERBOSE=0):
    '''Build a better consensus from an assembly of a few mapped reads
    
    This method exploits local linkage information to get frameshifts right.
    '''
    from hivwholeseq.utils.mapping import extract_mapped_reads_subsample_object

    if VERBOSE >= 1:
        print 'Building consensus from mapped reads:', bamfilename

    # Collect a few reads from every block_len nucleotides
    reads_grouped = [[] for i in xrange(len_fragment / block_len)]
    poss = np.array([i * block_len for i in xrange(len_fragment / block_len)], int)

    read_pairs = extract_mapped_reads_subsample_object(bamfilename,
                                                       n_reads=maxreads,
                                                       maxreads=max(maxreads, 100000),
                                                       VERBOSE=VERBOSE)

    if VERBOSE >= 2:
        print 'Grouping subsample of read pairs'

    # Classify according to START position
    for read_pair in read_pairs:
        for read in read_pair:
            ind = (read.pos <= poss).nonzero()[0]
            if len(ind) and (np.abs(read.isize) > \
                             poss[ind[0]] - read.pos + block_len + 20):

                # Forget about exact mapping position, CIGAR, etc. We redo this
                # via MSA. It is only required that the start position is covered
                reads_grouped[ind[0]].append(read.seq)
    
    if VERBOSE >= 4:
        print 'Reads per group:', map(len, reads_grouped)


    # Trim the last groups which have no reads
    # NOTE: even tiny groups are useful: of course, more errors there
    for i in xrange(len(reads_grouped) - 1, -1, -1):
        if len(reads_grouped[i]) < 2:
            del reads_grouped[i]
    poss = poss[:len(reads_grouped)]

    if VERBOSE >= 2:
        print 'Aligning groups and picking consensus blocks'

    # MSA within each group
    conss = []
    for irg, read_grouped in enumerate(reads_grouped):

        # Pick a few random ones out (or all if there are less than that)
        ind = np.arange(len(read_grouped))
        np.random.shuffle(ind)
        ind = ind[:min_reads_per_group]
        read_grouped_rnd = [read_grouped[i] for i in ind]
        ali = align_muscle(*[SeqRecord(Seq(s, ambiguous_dna), id=str(i))
                             for (i, s) in enumerate(read_grouped_rnd)])

        # Trim alignment to start from the first common position (this exists
        # because we binned reads according to their start!). In addition, every
        # alignment must cover at least a block length plus some overlap
        alim = np.array(ali)
        pos_start = ((alim != '-').mean(axis=0) > 0.9).nonzero()[0][0]
        # The end of the fragment should be a hard limit for all reads, because
        # of the primer trimming
        if irg == len(reads_grouped) - 1:
            pos_end = alim.shape[1]
        else:
            pos_end = pos_start + block_len + 20
        ali_trim = ali[:, pos_start: pos_end]

        if VERBOSE >= 4:
            print ali_trim
        
        # Get the most common string in this block (HERE WE USE LINKAGE INFO!)
        # NOTE: we could be more sophisticated here (cluster the reads, like in
        # a real-world assembly) but we should be needing this
        cou = Counter([str(s.seq) for s in ali])
        conss.append(str(Seq(cou.most_common(1)[0][0], ambiguous_dna).ungap('-')))    

    if VERBOSE >= 2:
        print 'Joining consensus blocks'
        
    # Join blocks (from the left)
    consensus = [conss[0]]
    for j, cons in enumerate(conss[1:], 1):
        seed = consensus[-1][-20:]
        sl = len(seed)
        pos_start = cons.find(seed)
        # Allow imperfect matches
        if pos_start == -1:
            consm = np.fromstring(cons, 'S1')
            seedm = np.fromstring(seed, 'S1')
            n_matches = [(consm[i: i + sl] == seedm).sum()
                         for i in xrange(len(cons) - len(seed))]
            pos_start = np.argmax(n_matches)

            # Try to only add non-bogus stuff
            if n_matches[pos_start] < 0.66 * sl:
                pos_start = -1
                if VERBOSE >= 4:
                    print 'block n.', j, 'not joint!'

        if pos_start != -1:
            consensus.append(cons[pos_start + sl:])
        
    consensus = ''.join(consensus)

    return consensus


def get_allele_frequencies_alignment(ali, alpha=alpha, VERBOSE=0):
    '''Get allele frequencies from MSA'''
    if len(ali.shape) > 1:
        af = np.zeros((len(alpha), len(ali[1])))
    else:
        af = np.zeros(len(alpha))
    alim = np.asarray(ali)
    for ia, a in enumerate(alpha):
        af[ia] = (alim == a).mean(axis=0)

    return af


def get_entropy(afs, alphabet_axis=None, VERBOSE=0):
    '''Get entropy from allele freqs'''
    if alphabet_axis is None:
        if len(afs.shape) == 1:
            alphabet_axis = 0
        else:
            alphabet_axis = -2

    S = np.maximum(0, -((afs) * np.log2(np.maximum(afs, 1e-8))).sum(axis=alphabet_axis))
    return S


# PLOT
def plot_coverage(data_folder, adaID, fragment, counts, VERBOSE=0, savefig=False):
    '''Plot figure with the coverage'''
    from hivwholeseq.sequencing.filenames import get_coverage_figure_filename as gff

    if VERBOSE >= 1:
        print 'Plotting coverage: '+adaID+' '+fragment

    coverage = counts.sum(axis=1).sum(axis=0)

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(15, 8))
    ax.plot(coverage + 0.5)
    ax.set_yscale('log')
    ax.set_xlabel('Position')
    ax.set_ylabel('Coverage')
    ax.set_title('adaID '+adaID+', fragment '+fragment)

    if savefig:
        outputfile = gff(data_folder, adaID, fragment)
        fig.savefig(outputfile)
        plt.close(fig)

    else:
        plt.ion()
        plt.show()



def plot_SFS_folded(data_folder, adaID, fragment, nu_filtered, VERBOSE=0, savefig=False):
    '''Plot the site frequency spectrum (folded)'''
    if VERBOSE >= 1:
        print 'Plotting folded SFS'

    from hivwholeseq.sequencing.filenames import get_SFS_figure_filename as gff
    import matplotlib.pyplot as plt
    import numpy as np

    nu_maj = np.ma.masked_all(nu_filtered.shape[1])
    nu_min = np.ma.masked_all(nu_filtered.shape[1])
    for pos, nus in enumerate(nu_filtered.T):
        if nus[0] == np.ma.masked:
            continue
        nus = np.sort(nus)
        if (nus[-1] < 0.5):
            if VERBOSE >= 3:
                print pos, 'has 3+ alleles:', nus, 'skipping.'
            continue

        nu_maj[pos] = nus[-1]
        nu_min[pos] = nus[-2]

    nu_maj_fold = 1 - nu_maj

    nu_mm = np.concatenate([nu_maj_fold, nu_min])
    nu_mm = np.array(nu_mm[nu_mm > 1e-5])
    nu_mm.sort()

    # Cumulative histogram
    fig, ax = plt.subplots(1, 1)
    ax.set_xlabel(r'$\nu$', fontsize=20)
    ax.set_ylabel('# alleles < x folded')
    ax.set_title('adaID '+adaID+', '+fragment)
    ax.set_xlim(10**(np.floor(np.log10(nu_mm[0] * 0.9))), 0.6)
    ax.set_xscale('log')
    ax.set_ylim(1.0 / len(nu_mm) * 0.9, 1.1)
    ax.set_yscale('log')
    ax.plot(nu_mm, 1.0 - np.linspace(0, 1 - 1.0 / len(nu_mm), len(nu_mm)), lw=2, c='b')

    if savefig:
        outputfile = gff(data_folder, adaID, fragment, cumulative=True, yscale='log')
        fig.savefig(outputfile)
        plt.close(fig)

    # Histogram
    fig, ax = plt.subplots(1, 1)
    ax.set_xlabel(r'$\nu$', fontsize=20)
    ax.set_ylabel('SFS folded (density)')
    ax.set_title('adaID '+adaID+', '+fragment)
    ax.set_xlim(10**(np.floor(np.log10(nu_mm[0] * 0.9))), 0.6)
    ax.set_xscale('log')
    ax.set_yscale('log')

    bins = np.logspace(-4, np.log10(0.5), 50)
    h = np.histogram(nu_mm, bins=bins, density=True)
    x = np.sqrt(h[1][1:] * h[1][:-1])
    y = h[0]
    ax.plot(x, y, lw=2, c='b')
    ax.scatter(x, y, s=50, edgecolor='none', facecolor='b')
    ax.grid()

    if savefig:
        outputfile = gff(data_folder, adaID, fragment, cumulative=False, yscale='log')
        fig.savefig(outputfile)
        plt.close(fig)
    else:
        plt.ion()
        plt.show()


