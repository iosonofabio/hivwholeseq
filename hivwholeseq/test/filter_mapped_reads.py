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

from miseq import alphal
from patients.filter_mapped_reads import filter_read_pair

from test.utils import Read


# Tests
class TestFilterReads(unittest.TestCase):

    def setUp(self):
        self.ref = 'CCCAAAGGGCCCTTTCCC'

        self.pairs = []

        # Trivial read pair, only matches
        readf = Read()
        readf.qname = 'matchonly'
        readf.pos = 3
        readf.seq = 'AAAGGGCCCTTT'
        readf.qual = 'G' * len(readf.seq)
        readf.cigar = [(0, len(readf.seq))]
        readf.is_reverse = False
        readf.is_unmapped = False
        readf.is_unpared = False
        readf.is_proper_pair = True

        readr = Read()
        readr.qname = 'matchonly'
        readr.pos = 5
        readr.seq = 'AGGGCCCTTTCCC'
        readr.qual = 'G' * len(readr.seq)
        readr.cigar = [(0, len(readr.seq))]
        readr.is_reverse = True
        readr.is_unmapped = False
        readr.is_unpared = False
        readr.is_proper_pair = True

        isize = readr.pos + len(readr.seq) - readf.pos
        readf.isize = isize
        readr.isize = -isize
        readf.mpos = readr.pos
        readr.mpos = readf.pos

        self.pairs.append((readf, readr))

        # The same, but reversed
        readf = Read()
        readf.qname = 'reverse'
        readf.pos = 3
        readf.seq = 'AAAGGGCCCTTT'
        readf.qual = 'G' * len(readf.seq)
        readf.cigar = [(0, len(readf.seq))]
        readf.is_reverse = False
        readf.is_unmapped = False
        readf.is_unpared = False
        readf.is_proper_pair = True

        readr = Read()
        readr.qname = 'reverse'
        readr.pos = 5
        readr.seq = 'AGGGCCCTTTCCC'
        readr.qual = 'G' * len(readr.seq)
        readr.cigar = [(0, len(readr.seq))]
        readr.is_reverse = True
        readr.is_unmapped = False
        readr.is_unpared = False
        readr.is_proper_pair = True

        isize = readr.pos + len(readr.seq) - readf.pos
        readf.isize = isize
        readr.isize = -isize
        readf.mpos = readr.pos
        readr.mpos = readf.pos

        self.pairs.append((readr, readf))
        
        # TODO: add more complex cases

        self.pairsdict = {reads[0].qname: reads for reads in self.pairs}


    def test_matchonly(self):
        '''Test match only pairs'''
        from copy import deepcopy

        pair = self.pairsdict['matchonly']
        pair_orig = deepcopy(pair)

        # Call the function
        pair_type = filter_read_pair(pair, self.ref, match_len_min=3)

        self.assertEqual(pair_type, 'good')
        self.assertEqual(pair[0], pair_orig[0])
        self.assertEqual(pair[1], pair_orig[1])






    #TODO: write a test, it has a bug if reads start before the fragment start



if __name__ == '__main__':
    unittest.main()

