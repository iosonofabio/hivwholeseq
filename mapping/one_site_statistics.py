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

from mapping.miseq import alpha, read_types
from mapping.mapping_utils import get_ind_good_cigars


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



