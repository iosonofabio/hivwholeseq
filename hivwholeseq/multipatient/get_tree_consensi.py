# vim: fdm=marker
'''
author:     Fabio Zanini
date:       22/06/14
content:    Take all consensi from sequenced samples connected to patient samples,
            make a multiple sequence alignment (MSA), and build a tree with it.

            Seq repetitions and PCR1/2 of the same sample should cluster, then within
            a patient. Else, contamination!
'''
# Modules
import os
import sys
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo

from hivwholeseq.patients.patients import load_patients, SamplePat
from hivwholeseq.patients.patients import load_samples_sequenced as lssp
from hivwholeseq.utils.mapping import align_muscle
from hivwholeseq.patients.filenames import get_consensi_alignment_filename, \
        get_consensi_tree_filename
from hivwholeseq.utils.tree import build_tree_fasttree
from hivwholeseq.reference import load_custom_reference
from hivwholeseq.utils.nehercook.ancestral import ancestral_sequences
from hivwholeseq.utils.tree import tree_to_json
from hivwholeseq.utils.generic import write_json



# Globals
refnames = {'38304': 'B',
            '38540': 'C',
            'LAI-III': 'B'}



# Functions
def find_feature_seq(record, region):
    '''Find feature in a sequence'''
    for feature in record.features:
        if feature.id == region:
            return feature.extract(record)
    else:
        raise ValueError('Region not found')


def trim_to_region(seq, regseq):
    '''Trim to a region via a similar sequence'''
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from seqanpy import align_overlap
    (score, ali1, ali2) = align_overlap(seq, regseq, score_gapopen=-20)
    start = len(ali2) - len(ali2.lstrip('-'))
    end = len(ali2.rstrip('-'))

    seqtrim = ali1[start: end].replace('-', '')

    return SeqRecord(Seq(seqtrim, seq.seq.alphabet),
                     id=seq.id,
                     name=seq.name,
                     description=seq.description)


def get_region_fragment_sample(sample, refreg, maxdelta=300, VERBOSE=0):
    '''Get the fragment of a region in an annotated reference sequence'''
    import numpy as np
    from seqanpy import align_overlap

    deltas = {}
    for fragment in ['F'+str(i) for i in xrange(1, 7)]:
        try:
            consensus_frag = sample.get_consensus(fragment)
        except IOError, OSError:
            if VERBOSE >= 2:
                print 'ERROR: consensus file not found'
            continue
        
        (score, ali1, ali2) = align_overlap(consensus_frag, refreg, score_gapopen=-20)
        start = len(ali2) - len(ali2.lstrip('-'))
        end = len(ali2.rstrip('-'))

        deltas[fragment] = (end - start) * 3 - score

    if len(deltas) == 0:
        return 'genomewide'

    fragments, deltas = zip(*deltas.iteritems())

    imax = np.argmin(deltas)
    if deltas[imax] > maxdelta:
        return 'genomewide'
    else:
        return fragments[imax]


def add_ancestral_subtype(tree, VERBOSE=0):
    '''Infer ancestral subtype from leaves'''
    from operator import attrgetter
    for node in tree.find_clades(terminal=False, order='postorder'):
        subtype_children = set(map(attrgetter('subtype'), node.clades))
        if (None in subtype_children) or (len(subtype_children) > 1):
            node.subtype = None
        else:
            node.subtype = subtype_children.pop()


def annotate_tree(tree, annotations, VERBOSE=0):
    '''Annotate tree'''
    from hivwholeseq.utils.tree import add_mutations_tree

    for leaf in tree.get_terminals():
        name = leaf.name
        
        if VERBOSE >= 2:
            print name

        if name not in annotations:
            continue

        anno = annotations[name]
        for key, value in anno.iteritems():
            setattr(leaf, key, value)

    add_mutations_tree(tree, translate=False,
                       sequence_attrname='sequence',
                       mutation_attrname='muts')

    add_ancestral_subtype(tree, VERBOSE=VERBOSE)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Build tree of all patient samples')
    pats_or_samples = parser.add_mutually_exclusive_group(required=False)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to map (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--regions', nargs='+',
                        help='Regions to analyze (e.g. V3 F6)')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--no-summary', action='store_false', dest='summary',
                        help='Do not save results in a summary file')
    #FIXME: there should be a --save argument

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    regions = args.regions
    submit = args.submit
    VERBOSE = args.verbose
    summary = args.summary


    fragments = ['F'+str(i) for i in xrange(1, 7)]

    patients = load_patients()

    refseq = load_custom_reference('HXB2', 'gb')

    # If no region specified, take all except tat and rev (because they are split)
    if regions is None:
        regions = [fea.id for fea in refseq.features
                   if len(fea.location.parts) == 1] + fragments

    # Get genomewide external references
    refs = []
    for refname, subtype in refnames.iteritems():
        ref = load_custom_reference(refname)
        ref.id = refname+'_'+subtype
        ref.name = refname+'_'+subtype
        ref.description = refname+' (subtype '+subtype+')'
        refs.append(ref)

    samples = lssp()
    # Collect all sequenced samples from patients
    if pnames is not None:
        samples = samples.loc[samples.patient.isin(pnames)]
        
    elif samplenames is not None:
        samples = samples.loc[samplenames]

    if VERBOSE >= 2:
        print 'samples', samples.index.tolist()

    for region in regions:
        consensi = []
        annotations = {}

        # Add external references
        if region in fragments:
            for refname, subtype in refnames.iteritems():
                ref = load_custom_reference(refname+'_'+region)
                ref.id = refname+'_'+refnames[refname]
                ref.name = refname+'_'+refnames[refname]
                ref.description = refname+' (subtype '+subtype+')'
                annotations[ref.id] = {'subtype': subtype}
                consensi.append(ref)

        else:
            refseqreg = find_feature_seq(refseq, region)
            for ref in refs:
                refreg = trim_to_region(ref, refseqreg)
                annotations[ref.id] = {'subtype': refnames[ref.id.split('_')[0]]}
                consensi.append(refreg)

        # Add patient samples
        for samplename, sample in samples.iterrows():
            if VERBOSE >= 1:
                print samplename, region,
                if VERBOSE == 1:
                    print ''

            sample = SamplePat(sample)
            
            if region in fragments:
                fragment = region
            else:
                fragment = get_region_fragment_sample(sample, refseqreg, VERBOSE=VERBOSE)

            try:

                # Get the whole fragment, then trim to region
                consensus_frag = sample.get_consensus(fragment)
            except IOError, OSError:
                if VERBOSE >= 2:
                    print 'ERROR: consensus file not found'
                continue

            if region in fragments:
                consensus = consensus_frag
            else:
                consensus = trim_to_region(consensus_frag, refseqreg)

            patcode = patients.loc[sample.patient, 'code']
            subtype = patients.loc[sample.patient, 'Subtype']
            consensus.id = patcode+'_'+str(sample['days since infection'])
            consensus.name = consensus.id
            consensi.append(consensus)

            annotations[consensus.id] = {'CD4': sample['CD4+ count'],
                                         'VL': sample['viral load'],
                                         'DSI': sample['days since infection'],
                                         'subtype': subtype,
                                         'patient': patcode,
                                        }

            if VERBOSE >= 2:
                print 'OK'

        if VERBOSE:
            print 'N consensi:', len(consensi)
                
        # Align
        if VERBOSE:
            print 'Align'
        ali = align_muscle(*consensi)
        ali_fn = get_consensi_alignment_filename('all', region)
        AlignIO.write(ali, ali_fn, 'fasta')

        # Make tree
        if VERBOSE:
            print 'Build tree'
        tree = build_tree_fasttree(ali_fn, VERBOSE=VERBOSE)


        if VERBOSE >= 2:
            print 'Infer ancestral sequences'
        a = ancestral_sequences(tree, ali, alphabet='ACGT-N', copy_tree=False,
                                attrname='sequence', seqtype='str')
        a.calc_ancestral_sequences()
        a.cleanup_tree()

        
        if VERBOSE >= 2:
            print 'Annotate leaves'
        annotate_tree(tree, annotations, VERBOSE=VERBOSE)

        if VERBOSE >= 2:
            print 'Ladderize tree'
        tree.ladderize()


        if VERBOSE >= 2:
            print 'Jsonify'
            tree_json = tree_to_json(tree.root,
                                     fields=('DSI', 'sequence', 'muts',
                                             'VL', 'CD4', 'confidence',
                                             'subtype', 'patient'),
                                    )

        if VERBOSE >= 2:
            print 'Save tree to file'
        tree_fn = get_consensi_tree_filename('all', region, format='json')
        write_json(tree_json, tree_fn, indent=1)
