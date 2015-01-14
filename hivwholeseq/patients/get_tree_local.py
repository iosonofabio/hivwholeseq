# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/12/14
content:    Check the V3 alignment in p4/15313, as it looks weird.
'''
# Modules
import os
import argparse
from operator import itemgetter, attrgetter
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from Bio import Phylo

from hivwholeseq.patients.patients import load_patient
from hivwholeseq.sequence_utils import align_muscle
from hivwholeseq.tree_utils import build_tree_fasttree
from hivwholeseq.argparse_utils import RoiAction
from hivwholeseq.patients.get_tree_consensi import annotate_tree



# Functions
def load_alignments(filename):
    '''Load alignments from website file'''
    import zipfile, zlib
    from Bio import AlignIO
    import StringIO

    alis = []
    with zipfile.ZipFile(filename, 'r') as zf:
        for fn in zf.namelist():
            f = StringIO.StringIO(zf.read(fn))
            ali = {'time': float(fn.split('_')[0]),
                   'ali': AlignIO.read(f, 'fasta')}
            alis.append(ali)

    return alis


def get_region_count_trajectories(patient, region, VERBOSE=0, countmin=5):
    '''Get haplotype trajectories in a region (from the website alignments)'''
    import numpy as np
    from hivwholeseq.website.filenames import get_precompiled_alignments_filename

    filename = get_precompiled_alignments_filename(patient.code, region)
    alis = load_alignments(filename)
    
    seqs_set = set()
    for ali in alis:
        seqs_set |= set([''.join(seq).replace('-', '')
                         for seq in ali['ali']
                         if int(seq.name.split('_')[1]) >= countmin])
    seqs_set = list(seqs_set)
    
    hct = np.zeros((len(seqs_set), len(alis)), int)
    for it, ali in enumerate(alis):
        for seq in ali['ali']:
            s = ''.join(seq).replace('-', '')
            count = int(seq.name.split('_')[1])
            if count < countmin:
                continue
            iseq = seqs_set.index(s)
            hct[iseq, it] = count

    seqs_set = np.array(seqs_set, 'S'+str(np.max(map(len, seqs_set))))
    times = np.array(map(itemgetter('time'), alis))
    ind = np.array([i for i, t in enumerate(patient.times) if t in times])

    # Filter out all time points without any counts
    ind_keep = hct.any(axis=0)
    ind = ind[ind_keep]
    hct = hct[:, ind_keep]

    return (hct.T, ind, seqs_set)


def expand_duplicates_annotate_tree(tree, htf, times, seqs, minfreq=0.01):
    '''Annotate tree with days and colors'''
    from operator import attrgetter
    from Bio.Phylo.BaseTree import Clade

    # Reroot
    iseqroot = htf[:, 0].argmax()
    for leaf in tree.get_terminals():
        iseq = int(leaf.name[3:]) - 1
        if iseq == iseqroot:
            break
    else:
        raise ValueError('Not able to reroot!')
    tree.root_with_outgroup(leaf)

    # Expand duplicates (same sequence but different time points)
    leaves = tree.get_terminals()
    for leaf in leaves:
        iseq = int(leaf.name[3:]) - 1
        seq = seqs[iseq]
        ht = htf[iseq]

        # Check how many copies of the leaf we need
        indseq = set((ht > minfreq).nonzero()[0])
        indseq.add(ht.argmax())
        indseq = sorted(indseq)

        # One copy only, just rename
        if len(indseq) == 1:
            it = indseq[0]
            h = ht[it]
            t = '{:1.1f}'.format(times[it])

            leaf.name = str(t)+'_'+seq
            leaf.frequency = h
            leaf.sequence = seq
            leaf.DSI = t

            continue

        # Several copy, transform into internal node and append subtree
        for it in indseq:
            it = indseq[0]
            h = ht[it]
            t = '{:1.1f}'.format(times[it])

            name = str(t)+'_'+seq
            newleaf = Clade(branch_length=1e-4, name=name)
            newleaf.frequency = h
            newleaf.sequence = seq
            newleaf.DSI = t

            leaf.clades.append(newleaf)

        leaf.name = None
        leaf.confidence = 1.0


def annotate_tree_for_plot(tree, htf, times, seqs, minfreq=0.02):
    '''Add annotations for plotting'''
    from matplotlib import cm
    cmap = cm.jet

    last_tp = max(leaf.DSI for leaf in tree.get_terminals())
    def get_color(node):
        return map(int, np.array(cmap(node.DSI/last_tp*0.9)[:-1]) * 255)

    # Annotate leaves
    for leaf in tree.get_terminals():
        leaf.color = get_color(leaf)

        if leaf.hf >= minfreq:
            leaf.label = 't = '+str(int(leaf.DSI))+', f = '+'{:1.2f}'.format(leaf.frequency)
        else:
            leaf.label = ''

    # Color internal branches
    for node in tree.get_nonterminals(order='postorder'):
        node.label = ''
        node.DSI = np.mean([c.DSI for c in node.clades])
        node.color = get_color(node)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get local trees',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--roi', required=True, action=RoiAction,
                        help='Region of interest (e.g. F1 300 350 or V3 0 +oo)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of reads analyzed per sample')
    parser.add_argument('--plot', action='store_true',
                        help='Plot local haplotype trajectories')
    parser.add_argument('--freqmin', type=int, default=0.01,
                        help='Minimal frequency to keep the haplotype')
    parser.add_argument('--save', action='store_true',
                        help='Save alignment to file')

    args = parser.parse_args()
    pname = args.patient
    roi = args.roi
    VERBOSE = args.verbose
    maxreads = args.maxreads
    use_plot = args.plot
    freqmin = args.freqmin
    use_save = args.save

    patient = load_patient(pname)

    if VERBOSE >= 1:
        print 'Get haplotypes'
    try:
        region = roi[0]
        (hct, ind, seqs) = get_region_count_trajectories(patient, region,
                                                         VERBOSE=VERBOSE)
    except IOError:
        region = ' '.join(map(str, roi))
        (hct, ind, seqs) = patient.get_local_haplotype_count_trajectories(*roi,
                                                                          VERBOSE=VERBOSE)
    
    hct = hct.T
    htf = 1.0 * hct / hct.sum(axis=0)

    # Get initial top haplo
    iseq0 = htf[:, 0].argmax()
    seq0 = seqs[iseq0]

    # Exclude too rare haplos
    indseq = (htf >= freqmin).any(axis=1)
    seqs = seqs[indseq]
    htf = htf[indseq]

    times = patient.times[ind]

    if VERBOSE >= 1:
        print 'Align sequences'
    seqsali = align_muscle(*seqs, sort=True)

    if VERBOSE >= 1:
        print 'Build local tree'
    tree = build_tree_fasttree(seqsali, VERBOSE=VERBOSE)
    expand_duplicates_annotate_tree(tree, htf, times, map(''.join, seqsali), minfreq=0.01)
    tree.ladderize()

    if use_save:
        if VERBOSE >= 2:
            print 'Save tree (JSON)'
        from hivwholeseq.tree_utils import tree_to_json
        from hivwholeseq.generic_utils import write_json
        annotate_tree(patient, tree, ali=[seqsali[iseq0]], VERBOSE=VERBOSE)
        tree_json = tree_to_json(tree.root)
        fn_out = patient.get_local_tree_filename(region, format='json')
        write_json(tree_json, fn_out)

    
    if use_plot:

        annotate_tree_for_plot(tree, htf, times, seqs, minfreq=0.1)

        if VERBOSE >= 1:
            print 'Plot'
        fig, axs = plt.subplots(1, 2, figsize=(10, 5))
        fig.suptitle(patient.code+', '+region)
 
        # Plot haplotype frequencies
        ax = axs[1]
        for hf in htf:
            ax.plot(times, hf, lw=2)
        
        ax.set_xlabel('Time [days]')
        ax.set_ylabel('Haplotype frequency')
        ax.grid(True)

        # Plot tree
        ax = axs[0]
        Phylo.draw(tree, axes=ax, do_show=False, label_func=attrgetter('label'),
                   show_confidence=False)
        ax.grid(True)
        ax.set_ylim(ax.get_ylim()[0] * 1.04, -ax.get_ylim()[0] * 0.04)
        
        plt.tight_layout(rect=(0, 0, 1, 0.95))
        plt.ion()
        plt.show()
