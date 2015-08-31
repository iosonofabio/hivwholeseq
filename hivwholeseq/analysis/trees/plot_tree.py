# vim: fdm=indent
'''
author:     Fabio Zanini
date:       09/04/15
content:    Build a tree of a local haplotype, af various kinds.
            - normal: use all sites
            - syn: using only "synonymous sites". It should give deeper terminal
            branches than usually, because of recombination.
'''
# Modules
import os
import sys
import argparse
from itertools import izip
from operator import itemgetter, attrgetter
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.utils.tree import build_tree_fasttree
from hivwholeseq.utils.argparse import PatientsAction
from hivwholeseq.store.store_tree_consensi import annotate_tree
from hivwholeseq.utils.nehercook.ancestral import ancestral_sequences
from hivwholeseq.utils.tree import tree_to_json
from hivwholeseq.utils.generic import write_json

from hivwholeseq.store.store_haplotypes_alignment_tree import (
    expand_annotate_alignment, annotate_tree_time_freq_count,
    annotate_tree_for_plot)


# Globals
pnames = ['p11']
regions = ['p15']


# Functions
def get_syn_sites(seq):
    '''Get position of "synonymous sites"'''
    from hivwholeseq.utils.sequence import translate_with_gaps
    prot = translate_with_gaps(seq)
    pos = np.zeros(len(seq), bool)
    for aa in ['A', 'G', 'I', 'L', 'P', 'S', 'R', 'T', 'V']:
        pos[(prot == aa).nonzero()[0] * 3 + 2] = True

    return pos



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get local trees',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', action=PatientsAction, default=pnames,
                        help='Patients to analyze')
    parser.add_argument('--regions', nargs='+', default=regions,
                        help='Genomic region (e.g. IN or V3)')
    parser.add_argument('--verbose', type=int, default=2,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot local haplotype trajectories')
    parser.add_argument('--freqmin', type=int, default=0.01,
                        help='Minimal frequency to keep the haplotype')
    parser.add_argument('--save', action='store_true',
                        help='Save alignment to file')
    parser.add_argument('--kinds', nargs='+', choices=['syn', 'normal'],
                        default=['normal'],
                        help='What kind of trees to plot')

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
            
            times = patient.times[ind]
            hct = hct.T
            hft = 1.0 * hct / hct.sum(axis=0)

            # Duplicate sequences for tree
            ali = expand_annotate_alignment(alim, hft, hct, times, freqmin=freqmin, VERBOSE=VERBOSE)

            # Exclude too rare haplos
            indseq = (hft >= freqmin).any(axis=1)
            alim = alim[indseq]
            hft = hft[indseq]

            for criterion in ['normal', 'syn']:

                if criterion == 'normal':
                    ali_new = ali
                elif criterion == 'syn':
                    # Keep only syn sites
                    cons0 = alim[hft[:, 0].argmax()]
                    try:
                        pos = get_syn_sites(cons0)
                    except ValueError:
                        continue
                    alim_new = alim[:, pos]
                    ali_new = []
                    for seq, seqm in izip(ali, alim_new):
                        seq = SeqRecord(Seq(''.join(seqm), seq.seq.alphabet),
                                        id=seq.id, name=seq.name, description=seq.description)
                        ali_new.append(seq)
                    ali_new = MultipleSeqAlignment(ali_new)

                if VERBOSE >= 2:
                    print 'Build local tree'
                tree = build_tree_fasttree(ali_new, VERBOSE=VERBOSE)

                if VERBOSE >= 2:
                    print 'Infer ancestral sequences'
                a = ancestral_sequences(tree, ali_new, alphabet='ACGT-N', copy_tree=False,
                                        attrname='sequence', seqtype='str')
                a.calc_ancestral_sequences()
                a.cleanup_tree()

                if VERBOSE >= 2:
                    print 'Annotate with time and frequency'
                annotate_tree_time_freq_count(tree, ali_new)

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

                if use_plot:
                    annotate_tree_for_plot(tree, minfreq=0.1)

                    if VERBOSE >= 2:
                        print 'Plot'
                    fig, ax = plt.subplots()
                    ax.set_title(patient.code+', '+region+', '+criterion)
                    
                    Phylo.draw(tree, axes=ax, do_show=False, label_func=attrgetter('label'),
                               show_confidence=False)
                    ax.grid(True)
                    ax.set_ylim(ax.get_ylim()[0] * 1.04, -ax.get_ylim()[0] * 0.04)
                    
                    plt.tight_layout()
                    plt.ion()
                    plt.show()
