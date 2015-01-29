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
        self.read = Read('AAAGGGTTTCCC', pos=1)


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

        # Equality test (they are ints)
        np.testing.assert_array_equal(counts, counts_check)



class Insert1(unittest.TestCase):
    def setUp(self):
        read = Read('AAAGGGTTTCCC', pos=1)
        read.cigar = [(0, 3), (1, 3), (0, 6)]
        self.read = read


    def test(self):
        '''Test allele counts from read with an in-frame insertion'''
        counts = np.zeros((len(alphaal), (len(self.read.seq) // 3) + 3), int)

        # Expected result
        from Bio.Data.CodonTable import standard_dna_table
        t = standard_dna_table.forward_table
        counts_check = counts.copy()
        for pos in xrange(len(self.read.seq) // 3):
            if 1 <= pos < 2:
                continue
            elif pos >= 2:
                posc = pos - 1
            else:
                posc = pos
            
            aa = alphaal.index(t[self.read.seq[3 * pos: 3 * (pos + 1)]])
            counts_check[aa, posc] += 1

        # Call the function
        get_allele_counts_aa_read(self.read, 1, 90, counts)

        # Equality test (they are ints)
        np.testing.assert_array_equal(counts, counts_check)


class Insert2(unittest.TestCase):
    def setUp(self):
        read = Read('AAAGGGTTTCCC', pos=1)
        read.cigar = [(0, 3), (1, 2), (0, 7)]
        self.read = read


    def test1(self):
        '''Test allele counts from read with an out-of-frame insertion'''
        counts = np.zeros((len(alphaal), (len(self.read.seq) // 3) + 3), int)

        # Expected result
        from Bio.Data.CodonTable import standard_dna_table
        t = standard_dna_table.forward_table
        counts_check = counts.copy()
        codons = ['AAA', 'GTT', 'TCC']
        for posc, codon in enumerate(codons):
            aa = alphaal.index(t[codon])
            counts_check[aa, posc] += 1

        # Call the function
        get_allele_counts_aa_read(self.read, 1, 90, counts)

        # Equality test (they are ints)
        np.testing.assert_array_equal(counts, counts_check)


    def test2(self):
        '''Test allele counts from read with an out-of-frame insertion
        
        In addition, the protein starts out-of-frame with the read.

          A A A G T T T C C C
        0 1 2 3 4 5 6 7 8 9 0
        \___/ \___/ \___/ \__
          0     1     2     3

                     TTC
        '''
        counts = np.zeros((len(alphaal), (len(self.read.seq) // 3) + 3), int)

        # Expected result
        from Bio.Data.CodonTable import standard_dna_table
        t = standard_dna_table.forward_table
        counts_check = counts.copy()
        codons = [(2, 'TTC')]
        for posc, codon in codons:
            aa = alphaal.index(t[codon])
            counts_check[aa, posc] += 1

        # Call the function
        get_allele_counts_aa_read(self.read, 0, 90, counts)

        # Equality test (they are ints)
        np.testing.assert_array_equal(counts, counts_check)


class Indel1(unittest.TestCase):
    def setUp(self):
        read = Read('AAAGGGTTTCCC', pos=1)
        read.cigar = [(0, 2), (1, 2), (2, 5), (0, 8)]
        self.read = read


    def test(self):
        '''Test allele counts from read with indels
        
          A A - - - - - G G T T T C C C
        0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6
        \___/ \___/ \___/ \___/ \___/ \__
          0     1     2     3     4     5

               ---         GTT   TCC
        
        '''
        counts = np.zeros((len(alphaal), (len(self.read.seq) // 3) + 3), int)

        # Expected result
        from Bio.Data.CodonTable import standard_dna_table
        t = standard_dna_table.forward_table
        counts_check = counts.copy()
        codons = [(1, '---'), (3, 'GTT'), (4, 'TCC')]
        for posc, codon in codons:
            if codon == '---':
                aa = alphaal.index('-')
            else:
                aa = alphaal.index(t[codon])
            counts_check[aa, posc] += 1

        # Call the function
        get_allele_counts_aa_read(self.read, 0, 90, counts)

        # Equality test (they are ints)
        np.testing.assert_array_equal(counts, counts_check)



if __name__ == '__main__':
    unittest.main()

