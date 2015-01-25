# vim: fdm=indent
'''
author:     Fabio Zanini
date:       23/01/15
content:    Test suite for allele counts calls from the reads.
'''
# Modules
# NOTE: in theory this is not necessary?
import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))


import unittest
from collections import defaultdict, Counter
import numpy as np
from copy import deepcopy

from miseq import alphal
from store.filter_mapped_reads import filter_read_pair

from test.utils import Read, fix_pair


# Tests
class TestFilterReads(unittest.TestCase):
    ref = 'CCCAAAGGGCCCTTTCCC'


class MatchOnly(TestFilterReads):
    '''Trivial read pair, only matches'''
    def setUp(self):
        readf = Read('AAAGGGCCCTTT', pos=3,
                     qname='matchonly')
        readr = Read('AGGGCCCTTTCCC', pos=5,
                     qname='matchonly',
                     is_reverse=True)
        fix_pair((readf, readr))
        self.pair = (readf, readr)

    def test(self):
        '''Test match only pairs'''
        pair_orig = deepcopy(self.pair)

        # Call the function
        pair_type = filter_read_pair(self.pair, self.ref, match_len_min=3)

        self.assertEqual(pair_type, 'good')
        self.assertEqual(self.pair, pair_orig)


class Reversed(TestFilterReads):
    '''Trivial but readf/r are reversed'''
    def setUp(self):
        readf = Read('AAAGGGCCCTTT', pos=3,
                     qname='matchonly')
        readr = Read('AGGGCCCTTTCCC', pos=5,
                     qname='matchonly',
                     is_reverse=True)
        fix_pair((readf, readr))
        self.pair = (readr, readf)

    def test(self):
        '''Test match only pairs'''
        pair_orig = deepcopy(self.pair)

        # Call the function
        pair_type = filter_read_pair(self.pair, self.ref, match_len_min=3)

        self.assertEqual(pair_type, 'good')
        self.assertEqual(self.pair, pair_orig)


class StartsBefore(TestFilterReads):
    '''Starts before 0'''
    def setUp(self):
        readf = Read('GGCCAACCCAAAGGG', pos=0,
                     qname='startsbefore')
        readf.cigar = [(1, 6), (0, len(readf.seq) - 6)]
        readr = Read('CCAAAGGGCCCTTT', pos=1,
                     qname='startsbefore',
                     is_reverse=True)
        self.pair = (readf, readr)


    def test(self):
        '''Test match only pairs'''
        pair = deepcopy(self.pair)
        pair[0].seq = pair[0].seq[7:]
        pair[0].qual = pair[0].qual[7:]
        pair[0].pos = 1
        pair[0].cigar = [(0, len(pair[0].seq))]
        pair[1].mpos = 1
        fix_pair(pair)

        # Call the function
        pair_type = filter_read_pair(self.pair, self.ref,
                                     trim_bad_cigars=1,
                                     match_len_min=3)

        self.assertEqual(pair_type, 'good')
        self.assertEqual(self.pair, pair)


# TODO: add more complex cases
# TODO: go back to the mapped reads and figure out what was wrong with those...



if __name__ == '__main__':
    unittest.main()

