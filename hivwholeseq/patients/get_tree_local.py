# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/12/14
content:    Get trees of haplotype alignments.
'''
# Modules
import os
import argparse
from operator import itemgetter, attrgetter
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from Bio import Phylo

from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.utils.sequence import align_muscle
from hivwholeseq.utils.tree import build_tree_fasttree
from hivwholeseq.utils.argparse import RoiAction
from hivwholeseq.store.store_tree_consensi import annotate_tree
from hivwholeseq.utils.nehercook.ancestral import ancestral_sequences
from hivwholeseq.utils.tree import tree_to_json, filter_rare_leaves
from hivwholeseq.utils.generic import write_json



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


def annotate_tree_for_plot(tree, minfreq=0.02):
    '''Add annotations for plotting'''
    from matplotlib import cm
    cmap = cm.jet

    last_tp = max(leaf.DSI for leaf in tree.get_terminals())
    def get_color(node):
        return map(int, np.array(cmap(node.DSI/last_tp*0.9)[:-1]) * 255)

    # Annotate leaves
    for leaf in tree.get_terminals():
        leaf.color = get_color(leaf)

        if leaf.frequency >= minfreq:
            leaf.label = ('t = '+str(int(leaf.DSI))+
                          ', f = '+'{:1.2f}'.format(leaf.frequency))
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
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--roi', required=True, action=RoiAction,
                        help='Region of interest (e.g. F1 300 350 or V3 0 +oo)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of reads analyzed per sample')
    parser.add_argument('--plot', action='store_true',
                        help='Plot local haplotype trajectories')
    parser.add_argument('--freqmin', type=float, default=0.01,
                        help='Minimal frequency to keep the haplotype')

    args = parser.parse_args()
    pnames = args.patients
    roi = args.roi
    VERBOSE = args.verbose
    maxreads = args.maxreads
    use_plot = args.plot
    freqmin = args.freqmin

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        if VERBOSE >= 1:
            print pname
    
        if os.path.isfile(patient.get_local_tree_filename(roi[0], format='json')):
            if VERBOSE >= 2:
                print 'Get tree'
            region = roi[0]
            tree = patient.get_local_tree(region)

        elif os.path.isfile(patient.get_local_tree_filename(' '.join(map(str, roi)), format='json')):
            if VERBOSE >= 2:
                print 'Get tree'
            region = ' '.join(map(str, roi))
            tree = patient.get_local_tree(region)

        else:
            raise IOError('Tree file not found')

        if VERBOSE >= 2:
            print 'Filter out too rare leaves'
        filter_rare_leaves(tree, freqmin, VERBOSE=VERBOSE)

        if use_plot:
            if VERBOSE >= 2:
                print 'Annotate tree for plotting'
            annotate_tree_for_plot(tree, minfreq=0.1)

            if VERBOSE >= 2:
                print 'Plot'
            fig, ax = plt.subplots()
            ax.set_title(patient.code+', '+region)
            
            Phylo.draw(tree, axes=ax, do_show=False, label_func=attrgetter('label'),
                       show_confidence=False)
            ax.grid(True)
            ax.set_ylim(ax.get_ylim()[0] * 1.04, -ax.get_ylim()[0] * 0.04)
            
            plt.tight_layout()
            plt.ion()
            plt.show()
