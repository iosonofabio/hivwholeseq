# vim: fdm=marker
'''
author:     Fabio Zanini
date:       24/10/13
content:    Build genome wide reference from the 6 fragments.
'''
# Modules
import os
import sys
import argparse
from itertools import izip
import numpy as np
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO as SeqIO

from hivwholeseq.sequencing.samples import SampleSeq, load_sequencing_run
from hivwholeseq.patients.filenames import get_initial_reference_filename
from hivwholeseq.sequencing.check_overlaps import get_overlapping_fragments, get_overlap



# Functions
def get_overlap(c1, c2, VERBOSE=0):
    '''Get overlap between two subsequent fragment sequences'''
    c1 = ''.join(c1)
    c2 = ''.join(c2)

    # Find beginning of c2 in c1
    seed = c2[:30]
    sl = len(seed)
    pos = c1.rfind(seed)
    if pos != -1:
        start = pos
    else:
        c1m = np.fromstring(c1, 'S1')
        seed = np.fromstring(seed, 'S1')
        n_match = [(seed == c1m[i: i+sl]).sum()
                   for i in xrange(len(c1) - sl)]
        pos = len(c1) - sl - 1 - np.argmax(n_match[::-1])
        if n_match[pos] < 0.66 * sl:
            raise ValueError('Overlap not found!')
        else:
            start = pos

    # Find end of c1 in c2
    end = len(c1) - start
    ov1 = c1[start:]
    ov2 = c2[:end]
    if VERBOSE >= 3:
        print ov1; print ov2

    ov1m = np.fromstring(ov1, 'S1')
    ov2m = np.fromstring(ov2, 'S1')
    n_diff = (ov1m != ov2m).sum()
    if n_diff > 15:
        raise ValueError('Fragments not compatible!')

    return (start, end)


def merge_fragments(pname, VERBOSE=0):
    '''Merge references at overlapping pairs'''
    import warnings

    references = {}
    for frag in ['F'+str(i+1) for i in xrange(6)]:
        cons_fn = get_initial_reference_filename(pname, frag)
        if os.path.isfile(cons_fn):
            references[frag] = SeqIO.read(cons_fn, 'fasta')

    fragments = references.keys()
    pairs = get_overlapping_fragments(fragments)

    overlaps = {}
    for (frag1, frag2) in pairs:
        try:
            overlap = get_overlap(references[frag1], references[frag2], VERBOSE=VERBOSE)
            overlaps[(frag1, frag2)] = overlap
        except ValueError:
            print 'WARNING: '+frag1+' and '+frag2+' have different references.'
            raise

    reference = []
    fragments = sorted(fragments)
    for i, frag in enumerate(fragments):
        # If the start is not an overlap, start a new reference and copy all
        if (i == 0) or (fragments[i-1], frag) not in overlaps:
            cons = [[frag], str(references[frag].seq)]
            reference.append(cons)

        # copy from the end of the overlap on
        else:
            cons = reference[-1]
            cons[0].append(frag)
            tmp = overlaps[(fragments[i-1], frag)]
            if tmp is not None:
                (_, end) = tmp
                cons[1] = cons[1]+str(references[frag][end:].seq)
            else:
                cons[1] = cons[1]+('N' * 10)+str(references[frag].seq)

    # Make SeqRecords out of references
    for i, (frags, cons) in enumerate(reference):
        rec = SeqRecord(Seq(cons, IUPAC.ambiguous_dna),
                        id=pname, name=pname,
                        description=pname+', initial reference, '+'-'.join(frags))
        reference[i] = (frags, rec)

    return reference



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Merge references and allele frequencies')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    
    args = parser.parse_args()
    pname = args.patient
    VERBOSE = args.verbose

    reference = merge_fragments(pname, VERBOSE=VERBOSE)
    for (frags, cons) in reference:
        if len(frags) < 6:
            fn = '-'.join(frags)
        else:
            fn = 'genomewide'
        
        output_filename = get_initial_reference_filename(pname, fn)
        SeqIO.write(cons, output_filename, 'fasta')
        if VERBOSE >= 1:
            print 'Reference written:', fn

