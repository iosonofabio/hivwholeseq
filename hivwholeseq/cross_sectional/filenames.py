# vim: fdm=marker
'''
author:     Fabio Zanini
date:       12/01/15
content:    Support module for filenames.
'''
# Modules
from hivwholeseq.sequencing.filenames import reference_folder
from hivwholeseq.filenames import table_folder



# Functions
def get_ctl_epitope_map_filename():
    '''Get filename of epitope map from LANL'''
    fn = 'ctl_summary.csv'
    return table_folder+fn


def get_subtype_reference_alignment_filename(region, subtype='B',
                                             refname='HXB2',
                                             type='nuc',
                                             VERBOSE=0):
    '''Get the filename of a subtype reference alignment'''
    tree_ali_foldername = reference_folder+'alignments/pairwise_to_'+refname+'/'
    fn = tree_ali_foldername+region+'.'+subtype+'.'+type+'.aligned.fasta'
    if VERBOSE >= 3:
        print 'Alignment file:', fn
    return fn


def get_subtype_reference_alignment_consensus_filename(region, subtype='B',
                                                       refname='HXB2',
                                                       type='nuc',
                                                       VERBOSE=0,
                                                       format='fasta'):
    '''Get the filename of a subtype reference alignment allele frequencies'''
    tree_ali_foldername = reference_folder+'alignments/pairwise_to_'+refname+'/'
    fn = tree_ali_foldername+region+'.'+subtype+'.'+type+'.aligned_consensus.'+format
    if VERBOSE >= 3:
        print 'Consensus file:', fn
    return fn


def get_subtype_reference_alignment_tree_filename(region, subtype='B',
                                                  refname='HXB2',
                                                  type='nuc',
                                                  VERBOSE=0,
                                                  format='newick'):
    '''Get the filename of a subtype reference alignment allele frequencies'''
    tree_ali_foldername = reference_folder+'alignments/pairwise_to_'+refname+'/'
    fn = tree_ali_foldername+region+'.'+subtype+'.'+type+'.aligned_tree.'+format
    if VERBOSE >= 3:
        print 'Tree file:', fn
    return fn


def get_subtype_reference_alignment_allele_frequencies_filename(region, subtype='B',
                                                                refname='HXB2',
                                                                type='nuc',
                                                                VERBOSE=0):
    '''Get the filename of a subtype reference alignment allele frequencies'''
    tree_ali_foldername = reference_folder+'alignments/pairwise_to_'+refname+'/'
    fn = tree_ali_foldername+region+'.'+subtype+'.'+type+'.aligned_afs.npy'
    if VERBOSE >= 3:
        print 'Afs file:', fn
    return fn


def get_subtype_reference_alignment_entropy_filename(region, subtype='B',
                                                     refname='HXB2',
                                                     type='nuc',
                                                     VERBOSE=0):
    '''Get the filename of a subtype reference alignment entropy'''
    tree_ali_foldername = reference_folder+'alignments/pairwise_to_'+refname+'/'
    fn = tree_ali_foldername+region+'.'+subtype+'.'+type+'.aligned_entropy.npy'
    if VERBOSE >= 3:
        print 'Entropy file:', fn
    return fn


def get_subtype_reference_alignment_entropy_syn_filename(region, subtype='B',
                                                     refname='HXB2',
                                                     type='nuc',
                                                     VERBOSE=0):
    '''Get the filename of a subtype reference alignment syn entropy'''
    tree_ali_foldername = reference_folder+'alignments/pairwise_to_'+refname+'/'
    fn = tree_ali_foldername+region+'.'+subtype+'.'+type+'.aligned_entropy_syn.pickle'
    if VERBOSE >= 3:
        print 'Entropy file:', fn
    return fn


def get_raw_LANL_sequences_filename(region):
    '''Get filename of raw sequences from LANL'''
    return reference_folder+'raw_LANL/'+region+'.fasta'
