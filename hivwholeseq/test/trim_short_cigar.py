# vim: fdm=indent
'''
author:     Fabio Zanini
date:       24/01/15
content:    Tests for the functions to trim short CIGARs from mapped reads.
'''
# Modules
# NOTE: in theory this is not necessary?
import os, sys
#sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))


import unittest
from collections import defaultdict, Counter
import numpy as np
from copy import deepcopy

from hivwholeseq.miseq import alphal
from hivwholeseq.utils.mapping import trim_short_cigars, trim_short_cigars_pair
from hivwholeseq.test.utils import Read, fix_pair



# Tests
class TestTrimCigars(unittest.TestCase):
    pass


class T1(TestTrimCigars):
    def setUp(self):
        read = Read(seq=''.join(chr(i) for i in xrange(65, 91)))
        read.cigar = [(1, 3), (0, 3), (2, 45), (1, 2), (0, 13), (1, 3)]
        read.cigar.append((0,
                           (len(read.seq) -
                            sum(bl for (bt, bl) in read.cigar if bt in (0, 1)))
                          ))

        self.read = read


    def test_nopad(self):
        read = Read(seq=self.read.seq[8:8+13],
                    pos=48)
        read.cigar = [(0, 13)]

        trim_short_cigars(self.read, match_len_min=5, trim_pad=0)
        self.assertEqual(read, self.read)


    def test_pad(self):
        read = Read(seq=self.read.seq[8+2:8+13-2],
                    pos=50)
        read.cigar = [(0, 9)]

        trim_short_cigars(self.read, match_len_min=5, trim_pad=2)
        self.assertEqual(read, self.read)



if __name__ == '__main__':
    unittest.main()

