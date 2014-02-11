# vim: fdm=marker
'''
author:     Fabio Zanini
date:       05/02/14
content:    Utility functions for mapped reads.
'''
# Functions
def test_read_integrity(read, VERBOSE=True):
    '''Test integrity of a mapped read'''
    if read.pos < 0:
        if VERBOSE:
            print 'Read starts before 0 ('+str(read.pos)+'):', read.qname
        return True

    if ((not read.is_reverse) and (read.isize < 0)) or \
       (read.is_reverse and (read.isize > 0)):
        if VERBOSE:
            print 'Read not integer (sign of isize):', read.qname
        return True
    
    if abs(read.isize) != sum(bl for (bt, bl) in read.cigar if bt in (0, 2)):
        if VERBOSE:
            print 'Read not integer (insert size <-> CIGAR):', read.qname
        return True

    if sum(bl for (bt, bl) in read.cigar if bt in (0, 1)) != read.rlen:
        if VERBOSE:
            print 'Read not integer (CIGAR <-> seq):', read.qname
        return True

    return False
