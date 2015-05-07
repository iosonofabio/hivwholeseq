# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/12/14
content:    Store phylogenetic tree of local haplotypes/consensi.
'''
# Modules
import os
import sys
import argparse
from operator import itemgetter, attrgetter
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio import AlignIO

from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.utils.argparse import PatientsAction
from hivwholeseq.utils.tree import build_tree_fasttree
from hivwholeseq.store.store_tree_consensi import annotate_tree
from hivwholeseq.utils.nehercook.ancestral import ancestral_sequences
from hivwholeseq.utils.tree import tree_to_json
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


def expand_annotate_alignment(alim, hft, hct, times, freqmin=0.01, VERBOSE=0):
    '''Expand alignment to duplicate leaves that stay, and rename seqs'''
    from itertools import izip
    from hivwholeseq.utils.sequence import convert_alim_to_biopython

    seqs = []
    names = []
    for i, seq in enumerate(alim):
        for time, hf, hc in izip(times, hft[i], hct[i]):
            if hf >= freqmin:
                name = str(time)+'_days_'+str(int(100 * hf))+'%_'+str(hc)+'_1'
                while name in names:
                    name = '_'.join(name.split('_')[:-1]+
                                    [str(int(name.split('_')[-1]) + 1)])

                seqs.append(seq)
                names.append(name)

    ali = convert_alim_to_biopython(seqs)
    for name, seq in izip(names, ali):
        seq.id = seq.name = name
        seq.description = ''
        
    return ali


def annotate_tree_time_freq_count(tree, ali):
    '''Annotate tree with days and colors'''
    from operator import attrgetter

    # Assign times and freqs to leaves
    for leaf in tree.get_terminals():
        leaf.DSI = float(leaf.name.split('_')[0])
        leaf.frequency = float(leaf.name.split('_')[2].split('%')[0]) / 100.0
        leaf.count = int(leaf.name.split('_')[3])

    
    # Reroot
    tmin = min(leaf.DSI for leaf in tree.get_terminals())
    newroot = max((leaf for leaf in tree.get_terminals() if leaf.DSI == tmin),
                   key=attrgetter('frequency'))
    tree.root_with_outgroup(newroot)


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


def extract_alignment(tree, VERBOSE=0):
    '''Extract aligned sequences from phylogenetic tree, including duplicates'''
    from Bio.Align import MultipleSeqAlignment as MSA
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet.IUPAC import ambiguous_dna

    ali = MSA([SeqRecord(Seq(leaf.sequence, ambiguous_dna),
                         id='days_'+str(int(leaf.DSI))+'_frequency_'+'{:2.0%}'.format(leaf.frequency),
                         name='days_'+str(int(leaf.DSI))+'_frequency_'+'{:2.0%}'.format(leaf.frequency),
                         description=('days since infection: '+str(int(leaf.DSI))+', '+
                                      'frequency: '+'{:2.0%}'.format(leaf.frequency)+', '+
                                      'n.reads: '+str(int(leaf.count))))
               for leaf in tree.get_terminals()])
    return ali


def plot_tree(tree, title=''):
    '''Plot the haplotype tree'''
    annotate_tree_for_plot(tree, minfreq=0.1)

    fig, ax = plt.subplots()
    
    Phylo.draw(tree,
               axes=ax,
               do_show=False,
               label_func=attrgetter('label'),
               show_confidence=False)
    ax.grid(True)
    ax.set_ylim(ax.get_ylim()[0] * 1.04, -ax.get_ylim()[0] * 0.04)
    if title:
        ax.set_title(title)
    
    plt.tight_layout()



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Get local trees',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', action=PatientsAction,
                        help='Patients to analyze')
    parser.add_argument('--regions', required=True, nargs='+',
                        help='Genomic region (e.g. IN or V3)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot local haplotype trajectories')
    parser.add_argument('--freqmin', type=int, default=0.01,
                        help='Minimal frequency to keep the haplotype')
    parser.add_argument('--save', action='store_true',
                        help='Save alignment to file')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    use_plot = args.plot
    freqmin = args.freqmin
    use_save = args.save

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    for region in regions:
        for pname, patient in patients.iterrows():
            patient = Patient(patient)

            if VERBOSE >= 1:
                print pname, region
        
            if VERBOSE >= 2:
                print 'Get haplotypes'

            (hct, ind, alim) = patient.get_haplotype_count_trajectory(region, aligned=True)
            
            if len(ind) == 0:
                if VERBOSE >= 2:
                    print 'Not time points found. Skip'
                continue

            times = patient.times[ind]
            hct = hct.T
            hft = 1.0 * hct / hct.sum(axis=0)

            # Duplicate sequences for tree
            ali = expand_annotate_alignment(alim, hft, hct, times, freqmin=freqmin, VERBOSE=VERBOSE)

            # Exclude too rare haplos
            indseq = (hft >= freqmin).any(axis=1)
            alim = alim[indseq]
            hft = hft[indseq]

            if VERBOSE >= 2:
                print 'Build local tree'
            tree = build_tree_fasttree(ali, VERBOSE=VERBOSE)

            if VERBOSE >= 2:
                print 'Infer ancestral sequences'
            a = ancestral_sequences(tree, ali, alphabet='ACGT-N', copy_tree=False,
                                    attrname='sequence', seqtype='str')
            a.calc_ancestral_sequences()
            a.cleanup_tree()

            if VERBOSE >= 2:
                print 'Annotate with time and frequency'
            annotate_tree_time_freq_count(tree, ali)

            # FIXME: for the amino acid mutations, we must make sure that we are
            # codon aligned (with codon_align). The problem is that sometimes
            # V3 has gaps of 1-2 nucleotides... at frequency 10%?!
            # NOTE: this might be due to compensating indels right outside V3,
            # as seen in cross-sectional alignments
            if VERBOSE >= 2:
                print 'Annotate tree'
            annotate_tree(patient, tree, VERBOSE=VERBOSE)
                
            if VERBOSE >= 2:
                print 'Ladderize tree'
            tree.ladderize()

            if use_save:
                if VERBOSE >= 2:
                    print 'Save tree (JSON)'
                fn = patient.get_local_tree_filename(region, format='json')
                tree_json = tree_to_json(tree.root,
                                         fields=('DSI', 'sequence',
                                                 'muts',
                                                 'VL', 'CD4',
                                                 'frequency',
                                                 'count',
                                                 'confidence'),
                                        )
                write_json(tree_json, fn)

            if VERBOSE >= 2:
                print 'Extract alignment from tree (so with duplicates)'
            ali_tree = extract_alignment(tree, VERBOSE=VERBOSE)

            if use_save:
                if VERBOSE >= 2:
                    print 'Save alignment from tree'
                fn = patient.get_haplotype_alignment_filename(region, format='fasta')
                AlignIO.write(ali_tree, fn, 'fasta')



    if use_plot:
        plt.ion()
        plt.show()
