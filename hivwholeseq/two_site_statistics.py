# vim: fdm=marker
'''
author:     Fabio Zanini
date:       20/03/14
content:    Collection of functions to do paired site statistics (coallele counts,
            cocoverage, correlations, linkage disequilibrium).
'''
# Modules
import sys
import numpy as np
from itertools import izip
import pysam

from hivwholeseq.miseq import alpha, alphal
from hivwholeseq.mapping_utils import pair_generator




# Functions
def get_coallele_counts_from_file(bamfilename, length, qual_min=35,
                                  maxreads=-1, VERBOSE=0,
                                  use_tests=False):
    '''Get counts of join occurence of two alleles'''
    # TODO: add quality thresholds (easy-peasy)

    if VERBOSE >= 1:
        print 'Getting coallele counts'
    
    if VERBOSE >= 2:
        print 'Initializing matrix of cocounts'
    # NOTE: ignore fwd/rev and read1/2 for now
    counts = np.zeros((len(alpha), len(alpha), length, length), int)
    posall = np.zeros(1000, dtype=[('pos', int), ('aind', int)])

    if VERBOSE >= 2:
        print 'Scanning reads'
    # NOTE: the reads should already be filtered of unmapped stuff at this point
    with pysam.Samfile(bamfilename, 'rb') as bamfile:

        for ir, reads in enumerate(pair_generator(bamfile)):
            if ir == maxreads:
                if VERBOSE:
                    print 'Max read number reached:', maxreads
                break
        
            if (VERBOSE >= 2) and (not ((ir +1) % 10)):
                if (VERBOSE == 2) and (ir + 1 != 10):
                    sys.stdout.write("\x1b[1A")
                print (ir+1) 

            # Check positions and CIGAR (this should REALLY not be necessary)
            if use_tests:
                for read in reads:
                    pos_ref = read.pos
                    pos_read = 0
                    if VERBOSE >= 3:
                        print read.pos, read.cigar, read.seq[:20]
                    for (bt, bl) in read.cigar:
                        if (bt in [0, 1, 2]) and (pos_ref > length):
                            raise ValueError('Pos exceeded the length of the fragment')
                    # Insertions
                    if bt == 1:
                        pos_read += bl
                    # Deletions
                    elif bt == 2:
                        pos_ref += bl
                    # Matches
                    elif bt == 0:
                        pos_read += bl
                        pos_ref += bl
                    # Other types of cigar?
                    else:
                        raise ValueError('CIGAR type '+str(bt)+' not recognized')

            # Temp structures
            posall[:] = (-1, -1)
            iall = 0

            # Collect alleles
            for read in reads:
                alleles_ind = np.fromiter((alphal.index(x) for x in read.seq), int, len(read.seq))
                pos_ref = read.pos
                pos_read = 0
                for (bt, bl) in read.cigar:
                    if bt == 1:
                        pos_read += bl
                    elif bt == 2:
                        posall['pos'][iall: iall + bl] = np.arange(pos_ref, pos_ref + bl)
                        posall['aind'][iall: iall + bl] = 4
                        iall += bl
                        pos_ref += bl
                    else:
                        alleles_indb = alleles_ind[pos_read: pos_read + bl]
                        for i in xrange(len(alpha)):
                            aitmp = (alleles_indb == i).nonzero()[0] + pos_ref
                            aitmplen = len(aitmp)
                            posall['pos'][iall: iall + aitmplen] = aitmp
                            posall['aind'][iall: iall + aitmplen] = i
                            iall += aitmplen

                        pos_read += bl
                        pos_ref += bl

            # Avoid doubles (paired reads are twice the same biological molecule)
            posall.sort(order=('pos', 'aind'))
            iall = (posall['pos'] != -1).nonzero()[0][0]
            while iall < len(posall) - 1:
                if posall['pos'][iall + 1] == posall['pos'][iall]:
                    # If the dups agree, skip one
                    if posall['aind'][iall + 1] == posall['aind'][iall]:
                        posall[iall + 1] = (-1, -1)

                    # else, pick one at random (FIXME: pick the highest phred)
                    else:
                        ibin = np.random.randint(2)
                        posall[iall + ibin] = (-1, -1)

                    iall += 1
                iall += 1

            # Add allele cocounts to the matrix
            # NOTE: this already takes care of the symmetry
            poss = [posall['pos'][posall['aind'] == i1] for i1 in xrange(len(alpha))]
            for i1 in xrange(len(alpha)):
                poss1 = poss[i1]
                if not len(poss1):
                    continue
                for i2 in xrange(len(alpha)):
                    poss2 = poss[i2]
                    if not len(poss2):
                        continue

                    # Raveling vodoo for efficiency - what python cannot, cobra can do!
                    cobra = counts[i1, i2].ravel()
                    ind = poss1.repeat(len(poss2)) * length + np.tile(poss2, len(poss1))
                    cobra[ind] += 1

    return counts
