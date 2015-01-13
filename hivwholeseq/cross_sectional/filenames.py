# vim: fdm=marker
'''
author:     Fabio Zanini
date:       12/01/15
content:    Support module for filenames.
'''
# Modules
from hivwholeseq.sequencing.filenames import reference_folder



# Functions
def get_subtype_reference_alignment_filename(region, subtype='B',
                                             refname='HXB2',
                                             type='nuc',
                                             VERBOSE=0):
    '''Get the filename of a subtype reference alignment'''
    tree_ali_foldername = reference_folder+'alignments/pairwise_to_'+refname+'/'
    ali_fn = tree_ali_foldername+region+'.'+subtype+'.'+type+'.aligned.fasta'
    if VERBOSE >= 3:
        print 'Alignment file:', ali_fn
    return ali_fn


def get_subtype_reference_alignment_entropy_filename(region, subtype='B',
                                                     refname='HXB2',
                                                     type='nuc',
                                                     VERBOSE=0):
    '''Get the filename of a subtype reference alignment'''
    tree_ali_foldername = reference_folder+'alignments/pairwise_to_'+refname+'/'
    ali_fn = tree_ali_foldername+region+'.'+subtype+'.'+type+'.aligned_entropy.npy'
    if VERBOSE >= 3:
        print 'Entropy file:', ali_fn
    return ali_fn


def get_subtype_reference_alignment_entropy_syn_filename(region, subtype='B',
                                                     refname='HXB2',
                                                     type='nuc',
                                                     VERBOSE=0):
    '''Get the filename of a subtype reference alignment'''
    tree_ali_foldername = reference_folder+'alignments/pairwise_to_'+refname+'/'
    ali_fn = tree_ali_foldername+region+'.'+subtype+'.'+type+'.aligned_entropy_syn.pickle'
    if VERBOSE >= 3:
        print 'Entropy file:', ali_fn
    return ali_fn



