# vim: fdm=marker
'''
author:     Fabio Zanini
date:       13/01/15
content:    Get the tree of consensi from a patient.
'''
# Modules
import sys
import os
import argparse
from operator import attrgetter
import numpy as np
from Bio import SeqIO, AlignIO

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.sequencing.samples import SampleSeq
from hivwholeseq.patients.patients import load_patients, Patient, SamplePat
from hivwholeseq.tree_utils import build_tree_fasttree
from hivwholeseq.utils.nehercook.ancestral import ancestral_sequences
from hivwholeseq.tree_utils import tree_to_json
from hivwholeseq.generic_utils import write_json



# Functions
def annotate_tree(patient, tree, VERBOSE=0,
                  fields=('DSI', 'muts', 'VL', 'ntemplates', 'CD4')):
    '''Annotate a tree with info on the nodes'''
    for node in tree.get_terminals():
        label = node.name
        entries = label.split('_')
        node.name = entries[0]

        if node.name == 'reference':
            continue

        # Days Since Infection
        if 'DSI' in fields:
            time = float(entries[0])
            node.DSI = time
        
        sample = patient.samples.loc[patient.times == time].iloc[0]
    
        if 'CD4' in fields:
            node.CD4 = sample['CD4+ count']

        if 'VL' in fields:
            node.VL = sample['viral load']

        # FIXME: shall we check the patient method for this?
        if 'ntemplates' in fields:
            node.ntemplates = sample['n templates']

    # For the mutations on branches, all nodes must have sequences already
    if ('muts') in fields or ('mutsprot' in fields):
        from Bio.Seq import translate
        from hivwholeseq.tree_utils import find_parent
        # NOTE: the root has no mutations by definition
        for node in tree.get_terminals() + tree.get_nonterminals():
            try:
                parent = find_parent(tree, node)
            except IndexError:
                continue

            pseq = parent.sequence
            nseq = node.sequence

            for attrname in ('muts', 'mutsprot'):
                if attrname not in fields:
                    continue
                
                if attrname == 'mutsprot':
                    pseq = translate(pseq)
                    nseq = translate(nseq)

                pseqm = np.fromstring(pseq, 'S1')
                nseqm = np.fromstring(nseq, 'S1')
                posm = (pseqm != nseqm).nonzero()[0]
                muts = [pseqm[p]+str(p+1)+nseqm[p] for p in posm]
                setattr(node, attrname, ', '.join(muts))



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Align consensi',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+', default=['all'],
                        help='Patients to analyze')
    parser.add_argument('--regions', nargs='*',
                        help='Regions to analyze (e.g. V3 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save alignment to file')
    parser.add_argument('--plot', action='store_true',
                        help='Plot phylogenetic tree. Requires --save and --tree')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    use_save = args.save
    use_plot = args.plot

    patients = load_patients()
    if pnames != ['all']:
        patients = patients.iloc[patients.index.isin(pnames)]

    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()

        if regions is None:
            refseq_gw = patient.get_reference('genomewide', 'gb')
            regionspat = map(attrgetter('id'), refseq_gw.features) + ['genomewide']
        else:
            regionspat = regions

        for region in regionspat:
            if VERBOSE >= 1:
                print pname, region
                if VERBOSE == 1:
                    print ''


            if VERBOSE >= 2:
                print 'Get alignment'
            ali = patient.get_consensi_alignment(region)


            if VERBOSE >= 2:
                print 'Build tree'
                sys.stdout.flush()
            tree = build_tree_fasttree(ali, rootname=ali[0].id,
                                       VERBOSE=VERBOSE)


            if VERBOSE >= 2:
                print 'Infer ancestral sequences'
            a = ancestral_sequences(tree, ali, alphabet='ACGT-N', copy_tree=False,
                                    attrname='sequence', seqtype='str')
            a.calc_ancestral_sequences()
            a.cleanup_tree()


            if VERBOSE >= 2:
                print 'Annotate tree'
            annotate_tree(patient, tree, VERBOSE=VERBOSE)


            if VERBOSE >= 2:
                print 'Ladderize tree'
            tree.ladderize()

            if use_save:
                if VERBOSE >= 2:
                    print 'Save tree (JSON)'
                fn = patient.get_consensi_tree_filename(region, format='json')
                tree_json = tree_to_json(tree.root,
                                         fields=('DSI', 'sequence', 'muts',
                                                 'VL', 'CD4', 'confidence'),
                                        )
                write_json(tree_json, fn)

            if use_plot:
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(figsize=(15, 12))
                Phylo.draw(tree, do_show=False, axes=ax)
                ax.set_title(pname+', '+region)

                x_max = max(tree.depths().itervalues())
                ax.set_xlim(0.995, 0.995 + (x_max - 0.995) * 1.4)
                ax.grid(True)
                
                plt.ion()
                plt.show()


