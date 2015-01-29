# vim: fdm=indent
'''
author:     Fabio Zanini
date:       29/01/15
content:    Test suite for amino acid allele counts calls from the reads.
'''
# Modules
# NOTE: in theory this is not necessary?
import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                os.pardir,
                                                os.pardir)))


import unittest
from collections import defaultdict, Counter
import numpy as np

from hivwholeseq.utils.sequence import alphaal
from hivwholeseq.one_site_statistics import get_allele_counts_aa_read

from hivwholeseq.test.utils import Read



# Tests
class MatchOnly(unittest.TestCase):
    def setUp(self):
        self.read = Read('AAAGGGTTTCCC', pos=1,
                         qname='matchonly')


    def test(self):
        '''Test allele counts from match only'''
        counts = np.zeros((len(alphaal), (len(self.read.seq) // 3) + 3), int)

        # Expected result
        from Bio.Data.CodonTable import standard_dna_table
        t = standard_dna_table.forward_table
        counts_check = counts.copy()
        for pos in xrange(len(self.read.seq) // 3):
            aa = alphaal.index(t[self.read.seq[3 * pos: 3 * (pos + 1)]])
            counts_check[aa, pos] += 1

        # Call the function
        get_allele_counts_aa_read(self.read, 1, 90, counts)

        print counts_check

        # Equality test (they are ints)
        np.testing.assert_array_equal(counts, counts_check)




if __name__ == '__main__':
    unittest.main()

