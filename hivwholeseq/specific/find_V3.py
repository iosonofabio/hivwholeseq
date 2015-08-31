# vim: fdm=marker
'''
author:     Fabio Zanini
date:       02/04/14
content:    Get reads from the V3 loop from any one sample for Thomas Leitner,
            with some data processing.
'''
# Modules
import argparse
import sys
import os
from itertools import izip
import numpy as np
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import ambiguous_dna

from seqanpy import align_global

from hivwholeseq.patients.patients import load_patient
from hivwholeseq.patients.filenames import get_mapped_to_initial_filename, \
        get_initial_reference_filename
from hivwholeseq.reference import load_HXB2
from hivwholeseq.utils.mapping import pair_generator




# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Extract seqs for Thomas')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    VERBOSE = args.verbose

    VERBOSE = 1
    pname = 'p1'
    itime = 0
    fragment = 'F5'
    maxgap = 0

    patient = load_patient(pname)

    cos = {}
    
    # Make sure it's a multiple of three, in case they translate
    # Plus the coordinates of Jan are a bit fuzzy, because the inner fwd primer
    # length is only 25, not 30...
    pos_V3_HXB2 = [6983, 7352]
    
    ref_rec = load_HXB2()
    refm = np.array(ref_rec)
    V3ref = ref_rec[pos_V3_HXB2[0]: pos_V3_HXB2[1]].seq
    
    cons_rec = patient.get_reference(fragment)
    cons = cons_rec.seq
    consm = np.array(cons_rec)
    
    ## Find the coordinates in the consensus
    V3primer_inner_fwd = np.fromstring('ACAATGYACACATGGAATTARGCCA', 'S1')
    seed = np.ma.array(V3primer_inner_fwd)
    seed[(seed == 'Y') | (seed == 'R')] = np.ma.masked
    sl = len(seed)
    n_matches = np.array([(consm[i: i + sl] == seed).sum() for i in xrange(len(consm) - sl)])
    start = (n_matches).argmax()
    if n_matches[start] < (sl - seed.mask.sum()) * 0.75:
        raise ValueError('Seed not found reliably')
    start += sl
    print 'V3 primer: fwd', start
    
    from Bio.Seq import reverse_complement as rc
    V3primer_inner_rev = np.fromstring(rc('AGAAAAATTCYCCTCYACAATTAAA'), 'S1')
    seed = np.ma.array(V3primer_inner_rev)
    seed[(seed == 'Y') | (seed == 'R')] = np.ma.masked
    sl = len(seed)
    n_matches = np.array([(consm[i: i + sl] == seed).sum() for i in xrange(len(consm) - sl)])
    end = (n_matches).argmax()
    if n_matches[end] < (sl - seed.mask.sum()) * 0.75:
        raise ValueError('Seed not found reliably')
    print 'V3 primer: rev', end
    
    V3con = cons[start: end]
    V3s = start
    V3e = end
    V3l = V3e - V3s
    print V3con.translate() 
