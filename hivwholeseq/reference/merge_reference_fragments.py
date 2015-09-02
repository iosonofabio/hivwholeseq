# vim: fdm=marker
'''
author:     Fabio Zanini
date:       16/01/15
content:    Join fragments for a reference (consensus) sequence into a
            genomewide sequence.
'''
# Modules
import os
import argparse

from hivwholeseq.reference import (
    load_custom_reference, get_custom_reference_filename)
from hivwholeseq.utils.sequence import merge_sequences



# Globals
fragments = ['F'+str(i) for i in xrange(1, 7)]



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Merge fragments in reference',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--reference', required=True,
                        help='Reference to analyze (e.g. LAI-III)')
    parser.add_argument('--fragments', nargs='+', default=fragments,
                        help='Fragments to merge')
    parser.add_argument('--save', action='store_true',
                        help='Save to file')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    refname = args.reference
    use_save = args.save
    VERBOSE = args.verbose


    consensi = [load_custom_reference(refname+'_'+fr) for fr in fragments]
    consensus = merge_sequences(consensi, VERBOSE=VERBOSE)

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    consrec = SeqRecord(Seq(consensus, alphabet=consensi[0].seq.alphabet),
                        id=refname+'_genomewide',
                        name=refname+'_genomewide',
                        description=refname+', genomewide reference (merged)',
                       )

    if use_save:
        if VERBOSE >= 1:
            print 'Save to file'
        from Bio import SeqIO
        fn = get_custom_reference_filename(refname, format='fasta')
        if os.path.isfile(fn):
            raise IOError('Destination file already exists')
        SeqIO.write(consrec, fn, 'fasta')

