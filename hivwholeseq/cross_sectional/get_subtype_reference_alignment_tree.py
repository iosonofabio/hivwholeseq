# vim: fdm=indent
'''
author:     Fabio Zanini
date:       18/03/15
content:    Make a tree of the reference alignment.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import ambiguous_dna, protein

from hivwholeseq.utils.sequence import alpha, alphaa
from hivwholeseq.utils.tree import build_tree_fasttree
from hivwholeseq.cross_sectional.filenames import (
    get_subtype_reference_alignment_tree_filename)
from hivwholeseq.cross_sectional.get_subtype_reference_alignment import (
    get_subtype_reference_alignment)



# Functions
def get_subtype_reference_alignment_tree(region,
                                         subtype='B',
                                         refname='HXB2',
                                         type='nuc',
                                         VERBOSE=0,
                                         format='newick'):
    '''Get filename of consensus of subtype reference alignment'''
    from Bio import Phylo
    from hivwholeseq.cross_sectional.filenames import (
        get_subtype_reference_alignment_tree_filename)

    fn = get_subtype_reference_alignment_tree_filename(region,
                                                       subtype=subtype,
                                                       refname=refname,
                                                       type=type,
                                                       VERBOSE=VERBOSE,
                                                       format=format)
    if format == 'newick':
        return Phylo.read(fn, format)
    
    from hivwholeseq.utils.tree import tree_from_json
    from hivwholeseq.utils.generic import read_json
    return tree_from_json(read_json(fn))




# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Calculate tree of subtype reference alignment',
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
                print 'Get alignment'
            ali = get_subtype_reference_alignment(region,
                                                  subtype=subtype,
                                                  refname=refname,
                                                  type=alitype,
                                                  VERBOSE=VERBOSE)

            if VERBOSE >= 2:
                print 'Calculate tree'
            tree = build_tree_fasttree(ali, VERBOSE=VERBOSE)


            if VERBOSE >= 2:
                print 'Save to file, newick'
            fn_out = get_subtype_reference_alignment_tree_filename(region,
                                                               subtype=subtype,
                                                               refname=refname,
                                                               type=alitype,
                                                               VERBOSE=VERBOSE,
                                                               format='newick')
            Phylo.write([tree], fn_out, 'newick')


            if VERBOSE >= 2:
                print 'Annotate tree (for JSON format)'
            from hivwholeseq.utils.nehercook.ancestral import ancestral_sequences
            a = ancestral_sequences(tree, ali, alphabet='ACGT-N', copy_tree=False,
                                    attrname='sequence', seqtype='str')
            a.calc_ancestral_sequences()
            del a


            if VERBOSE >= 2:
                print 'Save to file, JSON'
            fn_out = get_subtype_reference_alignment_tree_filename(region,
                                                               subtype=subtype,
                                                               refname=refname,
                                                               type=alitype,
                                                               VERBOSE=VERBOSE,
                                                               format='json')
            from hivwholeseq.utils.tree import tree_to_json
            from hivwholeseq.utils.generic import write_json
            tree_json = tree_to_json(tree.root, fields=('sequence', 'confidence'))
            write_json(tree_json, fn_out)


        else:
            consrec = get_subtype_reference_alignment_tree(region,
                                                           subtype=subtype,
                                                           refname=refname,
                                                           type=alitype,
                                                           VERBOSE=VERBOSE)
