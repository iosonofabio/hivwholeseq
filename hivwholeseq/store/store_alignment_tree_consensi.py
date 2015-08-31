# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/09/14
content:    Store alignment and phylogenetic tree of consensi.
'''
# Modules
import sys
import os
import argparse
from collections import defaultdict
from operator import attrgetter
import numpy as np
from Bio import SeqIO, AlignIO

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.reference import load_custom_reference
from hivwholeseq.patients.patients import load_patients, iterpatient, SamplePat
from hivwholeseq.utils.tree import build_tree_fasttree
from hivwholeseq.utils.mapping import align_muscle
from hivwholeseq.utils.nehercook.ancestral import ancestral_sequences
from hivwholeseq.utils.tree import tree_to_json, correct_minimal_branches
from hivwholeseq.utils.generic import write_json
from hivwholeseq.utils.argparse import PatientsAction
from hivwholeseq.genome_info import genes, proteins, RNA_structures
from hivwholeseq.patients.filenames import (
    get_consensi_alignment_filename,
    get_consensi_tree_filename)



# Globals
_regions = ['genomewide'] + list(genes) + list(proteins) + list(RNA_structures)
_refs = ['HXB2', 'NL4-3', '38304', '38540', 'F10']



# Functions
def annotate_tree(patient, tree, ali, VERBOSE=0):
    '''Annotate tree with metadata'''
    from hivwholeseq.utils.tree import add_mutations_tree
    from operator import attrgetter

    add_mutations_tree(tree, translate=False, mutation_attrname='muts')

    for leaf in tree.get_terminals():
        leaf.patient = patient.code
        leaf.DSI = float(leaf.name.split('_')[1])
        leaf.name = None

        sample = patient.samples.loc[patient.times == leaf.DSI].iloc[0]
        leaf.CD4 = sample['CD4+ count']
        leaf.VL = sample['viral load']
        leaf.subtype = patient.Subtype

    # Reroot
    tmin = min(leaf.DSI for leaf in tree.get_terminals())
    newroot = min((leaf for leaf in tree.get_terminals()),
                  key=attrgetter('DSI'))
    tree.root_with_outgroup(newroot)

    # Propagate the properties up the tree
    for prop in ['patient', 'subtype']:
        groups = defaultdict(list)
        for node in tree.get_nonterminals(order='postorder'):
            key = set([getattr(child, prop, None) for child in node.clades])
            if (None not in key) and (len(key) == 1):
                setattr(node, prop, key.pop())



def annotate_tree_all(tree, seqs, VERBOSE=0):
    '''Annotate tree with metadata'''
    from hivwholeseq.utils.tree import add_mutations_tree
    from operator import attrgetter

    add_mutations_tree(tree, translate=False, mutation_attrname='muts')

    seqnames = [seq.name for seq in seqs]

    for leaf in tree.get_terminals():
        if 'Isolate' not in leaf.name:
            leaf.patient = leaf.name.split('_')[0]
        else:
            leaf.patient = leaf.name[len('Isolate_'):]
        
        sample = seqs[seqnames.index(leaf.name)]
        leaf.subtype = sample.subtype

        if (leaf.patient[0] == 'p') and (len(leaf.patient) <= 3):
            leaf.DSI = float(leaf.name.split('_')[1])
            leaf.CD4 = sample.cell_count
            leaf.VL = sample.viral_load
        leaf.name = None

    # Reroot
    tree.root_at_midpoint()

    # Propagate the properties up the tree
    for prop in ['patient', 'subtype']:
        groups = defaultdict(list)
        for node in tree.get_nonterminals(order='postorder'):
            key = set([getattr(child, prop, None) for child in node.clades])
            if (None not in key) and (len(key) == 1):
                setattr(node, prop, key.pop())



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Align consensi and make trees',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', action=PatientsAction,
                        help='Patients to analyze')
    parser.add_argument('--regions', nargs='*', default=_regions,
                        help='Regions to analyze (e.g. V3 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save alignment to file')
    parser.add_argument('--joint', action='store_true',
                        help='Also make a joint alignment from selected patients')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    use_save = args.save
    use_joint = args.joint

    patients = load_patients()
    if pnames is not None:
        patients = patients.iloc[patients.index.isin(pnames)]

    pcodes = [p.code for _, p in patients.iterrows()]

    for region in regions:
        if use_joint:
            seqs_all = []
            for refname in _refs:
                seq = load_custom_reference(refname, region=region)
                if refname == 'F10':
                    refname = 'pZM246F_10'
                seq.id = 'Isolate_'+refname
                seq.name = 'Isolate_'+refname
                seq.description = 'Isolate: '+refname
                if refname in ['38540', 'pZM246F_10']:
                    seq.subtype = 'C'
                else:
                    seq.subtype = 'B'
                seqs_all.append(seq)

        for pname, patient in iterpatient(patients):
            patient.discard_nonsequenced_samples()
            if VERBOSE >= 1:
                print region, pname

            seqs = []
            for i, (sample) in enumerate(patient.itersamples()):
                if VERBOSE >= 2:
                    print sample.name,

                try:
                    cons_seq = sample.get_consensus(region)
                except IOError:
                    if VERBOSE >= 2:
                        print 'MISS'
                    continue

                time = str(int(sample['days since infection']))+'_days'
                cons_seq.id = patient.code+'_'+time
                cons_seq.name = cons_seq.id
                cons_seq.description = ', '.join(['Patient: '+patient.code,
                                                  time,
                                                  'region: '+region,
                                                  'consensus'])
                cons_seq.cell_count = sample['CD4+ count']
                cons_seq.viral_load = sample['viral load']
                cons_seq.subtype = patient['Subtype']
                seqs.append(cons_seq)
                if use_joint:
                    seqs_all.append(cons_seq)
                if VERBOSE >= 2:
                    print 'OK'
                    
            if VERBOSE >= 2:
                print 'Align',
            ali = align_muscle(*seqs, sort=True)
            if VERBOSE >= 2:
                print 'OK'

            if use_save:
                if VERBOSE >= 2:
                    print 'Save alignment',
                fn_out = patient.get_consensi_alignment_filename(region)
                mkdirs(os.path.dirname(fn_out))
                AlignIO.write(ali, fn_out, 'fasta')
                if VERBOSE >= 2:
                    print 'OK'

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
                print 'Annotate tree'
            annotate_tree(patient, tree, ali, VERBOSE=VERBOSE)

            if VERBOSE >= 2:
                print 'Ladderize tree'
            tree.ladderize()

            correct_minimal_branches(tree)

            if use_save:
                if VERBOSE >= 2:
                    print 'Save tree (JSON)'
                fn = patient.get_consensi_tree_filename(region, format='json')
                tree_json = tree_to_json(tree.root,
                                         fields=('DSI',
                                                 'patient',
                                                 'sequence',
                                                 'muts',
                                                 'VL',
                                                 'CD4',
                                                 'subtype',
                                                 'confidence'),
                                        )
                write_json(tree_json, fn, indent=1)

        
        if use_joint:
            if VERBOSE >= 2:
                print 'Align all patients',
            ali_all = align_muscle(*seqs_all, sort=True)
            if VERBOSE >= 2:
                print 'OK'

            if use_save:
                if VERBOSE >= 2:
                    print 'Save all patients',
                reg_tmp = '_'.join(pcodes)+'_'+region
                fn_out = get_consensi_alignment_filename('all', reg_tmp)
                mkdirs(os.path.dirname(fn_out))
                AlignIO.write(ali, fn_out, 'fasta')
                if VERBOSE >= 2:
                    print 'OK'

            if VERBOSE >= 2:
                print 'Build local tree'
            tree = build_tree_fasttree(ali_all, VERBOSE=VERBOSE)
            
            if VERBOSE >= 2:
                print 'Infer ancestral sequences'
            a = ancestral_sequences(tree, ali_all, alphabet='ACGT-N', copy_tree=False,
                                    attrname='sequence', seqtype='str')
            a.calc_ancestral_sequences()
            a.cleanup_tree()

            if VERBOSE >= 2:
                print 'Annotate tree'
            annotate_tree_all(tree, seqs_all, VERBOSE=VERBOSE)

            if VERBOSE >= 2:
                print 'Ladderize tree'
            tree.ladderize()

            correct_minimal_branches(tree)

            if use_save:
                if VERBOSE >= 2:
                    print 'Save tree (JSON)'
                reg_tmp = '_'.join(pcodes)+'_'+region
                fn = get_consensi_tree_filename('all', reg_tmp, format='json')
                tree_json = tree_to_json(tree.root,
                                         fields=('DSI',
                                                 'patient',
                                                 'sequence',
                                                 'muts',
                                                 'VL',
                                                 'CD4',
                                                 'subtype',
                                                 'confidence'),
                                        )
                write_json(tree_json, fn, indent=1)

