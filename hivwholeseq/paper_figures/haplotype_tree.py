# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/12/14
content:    Plot the quality along the read pair.
'''
# Modules
import os
import numpy as np

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.patients.patients import load_patient
from hivwholeseq.paper_figures.filenames import get_figure_folder
from hivwholeseq.sequence_utils import align_muscle
from hivwholeseq.tree_utils import build_tree_fasttree
from hivwholeseq.paper_figures.plots import plot_haplotype_tree_example



# Functions
def compress_data(tree, htf, pname, region):
    '''Compress data for plots, discarding useless info'''
    data = []
    datum = {'tree': tree,
             'htf': htf,
             'pname': pname,
             'region': region}
    data.append(datum)

    return data


def store_data(data, fn):
    '''Store data to file for the plots'''
    import cPickle as pickle
    with open(fn, 'wb') as f:
        pickle.dump(data, f, protocol=-1)


def load_data(fn):
    '''Load the data for the plots'''
    import cPickle as pickle
    with open(fn, 'rb') as f:
        return pickle.load(f)


def adjust_tree(tree, htf, times):
    '''Adjust tree for plotting, in place'''
    from Bio.Phylo.BaseTree import Clade

    # Reroot
    for leaf in tree.get_terminals():
        if int(leaf.name[3:]) - 1 == htf[0].argmax():
            tree.root_with_outgroup(leaf)
            break
    else:
        raise ValueError('root not found')

    # Duplicate leaves
    leaves = tree.get_terminals()
    for leaf in leaves:
        iseq = int(leaf.name[3:]) - 1
        ht = htf[:, iseq]

        # Check how many copies of the leaf we need
        indseq = np.sort(list(set((ht > 0.05).nonzero()[0]) | set([ht.argmax()])))

        # One copy only, just rename
        if len(indseq) == 1:
            it = indseq[0]
            h = ht[it]
            t = int(times[it])
            leaf.name = '{:1.0%}'.format(h)+' '+str(t)+' days'
            leaf.iseq = iseq
            leaf.it = it
            leaf.age = t
            leaf.age_frac = 1.0 * t / (times[-1] - times[0])
            continue

        # Several copy, transform into internal node and append subtree
        for it in indseq:
            it = indseq[0]
            h = ht[it]
            t = int(times[it])
            name = '{:1.0%}'.format(h)+' '+str(t)+' days'
            newleaf = Clade(branch_length=1e-4, name=name,
                            confidence=1.0)
            newleaf.age = t
            newleaf.age_frac = 1.0 * t / (times[-1] - times[0])
            newleaf.iseq = iseq
            newleaf.it = it
            leaf.clades.append(newleaf)

        leaf.name = None

    tree.ladderize()
    assign_age(tree)


def assign_age(tree):
    '''Assign an age to internal nodes recursively'''
    # NOTE: this should be onyl the root (which has a single child)
    if not hasattr(tree, 'clades'):
        tree.age_frac = assign_age(tree.clade)

    elif not hasattr(tree, 'age_frac'):
        tree.age_frac = np.mean([assign_age(node) for node in tree.clades])

    return tree.age_frac



# Script
if __name__ == '__main__':

    VERBOSE = 2

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    mkdirs(fn_data)
    fn_data = fn_data + 'haplotype_tree_example.pickle'

    if not os.path.isfile(fn_data):
        pname = '9669'
        region = 'V3'
        cutoff = 0.04

        patient = load_patient(pname)

        htc, ind, seqs = patient.get_region_count_trajectories(region)
        times = patient.times[ind]

        # Filter only seqs above the cutoff
        htf = (1.0 * htc.T / htc.sum(axis=1)).T
        indseq = (htf >= cutoff).any(axis=0)
        htc = htc[:, indseq]
        htf = htf[:, indseq]
        seqs = seqs[indseq]

        # Align, make tree
        ali = align_muscle(*seqs, sort=True)
        tree = build_tree_fasttree(ali, VERBOSE=VERBOSE)

        # Duplicate leaves and annotate tree
        adjust_tree(tree, htf, times)

        data = compress_data(tree, htf, pname, region)
        store_data(data, fn_data)
    else:
        data = load_data(fn_data)
        
    filename = foldername+'haplotype_tree_example'
    for ext in ['png', 'pdf', 'svg']:
        plot_haplotype_tree_example(data,
                                    VERBOSE=VERBOSE,
                                    savefig=filename+'.'+ext)

    plot_haplotype_tree_example(data,
                                VERBOSE=VERBOSE)
