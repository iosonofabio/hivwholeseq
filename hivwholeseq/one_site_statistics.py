# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/11/13
content:    Collection of functions to do single site statistics (allele counts,
            coverage, allele frequencies).
'''
# Modules
from collections import defaultdict
import numpy as np
import pysam

from hivwholeseq.miseq import alpha, read_types
from hivwholeseq.mapping_utils import get_ind_good_cigars


# Functions
def get_allele_counts_insertions_from_file(bamfilename, length, qual_min=35,
                                           maxreads=-1, VERBOSE=0):
    '''Get the allele counts and insertions'''
    # Prepare output structures
    counts = np.zeros((len(read_types), len(alpha), length), int)
    # Note: the data structure for inserts is a nested dict with:
    # position --> string --> read type --> count
    #  (dict)      (dict)       (list)      (int)
    inserts = defaultdict(lambda: defaultdict(lambda: np.zeros(len(read_types), int)))

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
        
            # Read CIGARs (they should be clean by now)
            seq = np.fromstring(read.seq, 'S1')
            qual = np.fromstring(read.qual, np.int8) - 33
            pos = read.pos

            # Iterate over CIGARs
            for ic, (block_type, block_len) in enumerate(read.cigar):

                # Check for pos: it should never exceed the length of the fragment
                if (block_type in [0, 1, 2]) and (pos >= length):
                    raise ValueError('Pos exceeded the length of the fragment')
            
                # Inline block
                if block_type == 0:
                    seqb = seq[:block_len]
                    qualb = qual[:block_len]
                    # Increment counts
                    for j, a in enumerate(alpha):
                        posa = ((seqb == a) & (qualb >= qual_min)).nonzero()[0]
                        if len(posa):
                            counts[js, j, pos + posa] += 1
            
                    # Chop off this block
                    if ic != len(read.cigar) - 1:
                        seq = seq[block_len:]
                        qual = qual[block_len:]
                        pos += block_len
            
                # Deletion
                elif block_type == 2:
                    # Increment gap counts
                    counts[js, 4, pos:pos + block_len] += 1
            
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
                        inserts[pos][seqb.tostring()][js] += 1
            
                    # Chop off seq, but not pos
                    if ic != len(read.cigar) - 1:
                        seq = seq[block_len:]
                        qual = qual[block_len:]
            
                # Other types of cigar?
                else:
                    raise ValueError('CIGAR type '+str(block_type)+' not recognized')

    return (counts, inserts)


def get_allele_counts_insertions_from_file_unfiltered(bamfilename, length, qual_min=30,
                                                      match_len_min=10,
                                                      maxreads=-1, VERBOSE=0):
    '''Get the allele counts and insertions'''
    # Prepare output structures
    counts = np.zeros((len(read_types), len(alpha), length), int)
    # Note: the data structure for inserts is a nested dict with:
    # position --> string --> read type --> count
    #  (dict)      (dict)       (list)      (int)
    inserts = defaultdict(lambda: defaultdict(lambda: np.zeros(len(read_types), int)))

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


def filter_nus(counts, coverage, VERBOSE=0):
    '''Filter allele frequencies from the four read types'''
    from scipy.stats import chi2_contingency

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

    # Set filtered allele frequencies
    nu_filtered = np.ma.masked_all((len(alpha), counts.shape[-1]))

    # Iterate over positions
    for i in xrange(counts.shape[-1]):
        
        # 1. if we cover neither fwd nor rev, keep masked
        if ind_low_cov_f[i] and ind_low_cov_b[i]:
            if VERBOSE >= 4:
                print 'pos', i, 'base', alpha[j], 'not covered'
            pass

        # 2. if we cover only one of them, well, just take the
        # arithmetic sum of counts
        elif ind_low_cov_f[i] != ind_low_cov_b[i]:
            nu_filtered[:, i] = 1.0 * counts[:, :, i].sum(axis=0) / coverage[:, i].sum()
            if VERBOSE >= 4:
                print 'pos', i, 'base', alpha[j], 'covered only once'
                
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
