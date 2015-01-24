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
from one_site_statistics import get_allele_counts_read

from test.utils import Read



# Tests
class TestAlleleCounts(unittest.TestCase):

    def setUp(self):
        self.reads = []

        # Trivial read, only matches
        read = Read('AAAGGGTTTCCC', pos=1,
                    qname='matchonly')
        self.reads.append(read)

        # Read that starts with an insertion
        read = Read('AAAGGGTTTCCC', pos=2,
                    qname='startins')
        read.cigar = [(1, 3), (0, len(read.seq) - 3)]
        self.reads.append(read)

        # TODO: add more complex cases

        self.readsdict = {read.qname: read for read in self.reads}


    def test_matchonly(self):
        '''Test allele counts from match only'''
        read = self.readsdict['matchonly']
        counts = np.zeros((len(alphal), len(read.seq) + read.pos + 1), int)
        inserts = defaultdict(lambda: Counter())

        # Expected result
        counts_check = counts.copy()
        for pos, nuc in enumerate(read.seq):
            counts_check[alphal.index(nuc), pos + read.pos] += 1

        # Call the function
        get_allele_counts_read(read, counts, inserts)

        # Equality test (they are ints)
        np.testing.assert_array_equal(counts, counts_check)


    def test_startins(self):
        '''Test allele counts from match only'''
        read = self.readsdict['startins']
        counts = np.zeros((len(alphal), len(read.seq) + read.pos + 1), int)
        inserts = defaultdict(lambda: Counter())

        # Expected result
        counts_check = counts.copy()
        for pos, nuc in enumerate(read.seq[read.cigar[0][1]:]):
            counts_check[alphal.index(nuc), pos + read.pos] += 1

        # Call the function
        get_allele_counts_read(read, counts, inserts)

        # Equality test (they are ints)
        np.testing.assert_array_equal(counts, counts_check)


#TODO: write a test class for the filtering function in /patients, because it
# has a bug with reads that start before the fragment start



if __name__ == '__main__':
    unittest.main()

