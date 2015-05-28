# vim: fdm=indent
'''
author:     Fabio Zanini
date:       04/02/15
content:    Calculate the subtype consensus.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import ambiguous_dna, protein

from hivwholeseq.utils.sequence import alpha, alphaa
from hivwholeseq.cross_sectional.filenames import (
    get_subtype_reference_alignment_allele_frequencies_filename,
    get_subtype_reference_alignment_consensus_filename)


# Functions
def get_subtype_reference_alignment_consensus(region,
                                              subtype='B',
                                              refname='HXB2',
                                              type='nuc',
                                              VERBOSE=0):
    '''Get filename of consensus of subtype reference alignment'''
    from Bio import SeqIO
    from hivwholeseq.cross_sectional.filenames import (
        get_subtype_reference_alignment_consensus_filename)

    fn = get_subtype_reference_alignment_consensus_filename(region,
                                                            subtype=subtype,
                                                            refname=refname,
                                                            type=type,
                                                            VERBOSE=VERBOSE)
    return SeqIO.read(fn, 'fasta')




# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Calculate consensus of subtype reference alignment',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Regions to analyze (e.g. F1 p17)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--subtype', default='B',
                        help='Subtype of the alignment')
    parser.add_argument('--reference', default='HXB2',
                        help='Reference of the alignment')
    parser.add_argument('--type', default='nuc',
                        help='nuc/aa nucleic or amino acid seqs')
    parser.add_argument('--save', action='store_true',
                        help='Save to file')

    args = parser.parse_args()
    regions = args.regions
    VERBOSE = args.verbose
    subtype = args.subtype
    refname = args.reference
    alitype = args.type
    use_save = args.save

    if alitype == 'nuc':
        alphabet = alpha
        alphabet_bio = ambiguous_dna
    else:
        alphabet = alphaa
        alphabet_bio = protein

    for region in regions:
        if VERBOSE >= 1:
            print region


        if use_save:
            if VERBOSE >= 2:
                print 'Get allele frequencies'
            fn_in = get_subtype_reference_alignment_allele_frequencies_filename(region,
                                                                      subtype=subtype,
                                                                      refname=refname,
                                                                      type=alitype,
                                                                      VERBOSE=VERBOSE)
            afs = np.load(fn_in)


            if VERBOSE >= 2:
                print 'Calculate consensus'
            
            consm = alphabet[afs[:5].argmax(axis=0)]
            consrec = SeqRecord(Seq(''.join(consm), alphabet_bio),
                                id='consensus_subtype'+subtype+'_refto_'+refname,
                                name='consensus_subtype'+subtype+'_refto_'+refname,
                                description=('Consensus of subtype '+subtype+
                                             ' reference alignment to '+refname),
                               )

            if VERBOSE >= 2:
                print 'Save to file'
            fn_out = get_subtype_reference_alignment_consensus_filename(region,
                                                               subtype=subtype,
                                                               refname=refname,
                                                               type=alitype,
                                                               VERBOSE=VERBOSE)
            SeqIO.write(consrec, fn_out, 'fasta')

        else:
            consrec = get_subtype_reference_alignment_consensus(region,
                                                               subtype=subtype,
                                                               refname=refname,
                                                               type=alitype,
                                                               VERBOSE=VERBOSE)
