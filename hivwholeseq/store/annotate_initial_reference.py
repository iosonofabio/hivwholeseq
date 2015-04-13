# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/09/14
content:    Annotate initial reference with positions of fragments, genes, other
            biologically relevant entities.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Alphabet.IUPAC import ambiguous_dna

from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.utils.sequence import correct_genbank_features_save




# Functions
def get_edges_fragments(patient, VERBOSE=0):
    '''Get the edges of fragments for this patient'''
    edges = {}
    for fragment in ['F'+str(i+1) for i in xrange(6)]:
        fn = patient.get_reference_filename(fragment)
        refseq = SeqIO.read(fn, 'fasta')
        if fragment == 'F1':
            edge_start = None
        else:
            edge_start = ''.join(refseq[:40])

        if fragment == 'F6':
            edge_end = None
        else:
            edge_end = ''.join(refseq[-40:])

        edge = (edge_start, edge_end)
        
        if VERBOSE >= 2:
            print fragment, str(edge_start), str(edge_end)

        edges[fragment] = edge

    return edges


def annotate_sequence(seqrecord, additional_edges={}, additional_features=['chunk'], VERBOSE=0):
    '''Annotate a consensus with the genes and stuff (in place)'''
    # TODO: what do we do with genes that do not start/end where they are
    # supposed to? Do we follow biology and track their new locations?
    from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
    from hivwholeseq.genome_info import gene_edges, RNA_structure_edges, \
            other_edges, find_region_edges, find_region_edges_multiple, \
            locate_gene
    edge_dict = {'gene': gene_edges,
                 'RNA structure': RNA_structure_edges,
                 'other': other_edges}
    edge_dict.update(additional_edges)
    additional_features = ['protein'] + additional_features
    features = edge_dict.keys() + additional_features

    if VERBOSE:
        print 'Features:', ', '.join(features)

    smat = np.array(seqrecord)

    for feature_type in edge_dict:
        edges_all = edge_dict[feature_type]
        print feature_type, edge_dict[feature_type].keys()
        for name, edges in edges_all.iteritems():
            if VERBOSE >= 2:
                print name,

            # Skip a feature if it's present already
            if name in map(lambda x: x.id, seqrecord.features):
                if VERBOSE >= 2:
                    print 'already present.'
                continue

            # Behave differently for unsplit regions and split ones
            if len(edges) == 2:
                # LTR problems with F6
                if 'F6' in name:
                    pos_edge = find_region_edges(smat[6000::], [edges[0], None])
                    pos_edge[0] += 6000
                elif feature_type == 'genes':
                    pos_edge = locate_gene(smat, name, output_compact=True)
                else:
                    pos_edge = find_region_edges(smat, edges)

                # Cut the primers for some features
                if (None not in pos_edge) and name in ['V1', 'V3', 'V4', 'V5']:
                    pos_edge[0] += len(edges[0])
                    pos_edge[1] -= len(edges[1])

                # Cut only the right primer for V2
                if (None not in pos_edge) and name in ['V2']:
                    pos_edge[1] -= len(edges[1])

                if pos_edge[0] is None:
                    if name not in ['F1', "LTR5'"]:
                        print 'WARNING: start not found'
                    pos_edge[0] = 0

                if pos_edge[1] is None:
                    if name not in ['F6', "LTR3'"]:
                        print 'WARNING: end not found'
                    pos_edge[1] = len(smat)

                location = FeatureLocation(*pos_edge)
            else:
                if feature_type == 'genes':
                    pos_edges = [locate_gene(smat, name+suff, output_compact=True)
                                 for suff in ('1', '2')]
                else:
                    pos_edges = find_region_edges_multiple(smat, edges)
                locations = [FeatureLocation(*pos_edge) for pos_edge in pos_edges]
                location = CompoundLocation(locations)

            if VERBOSE >= 2:
                print 'found:', location

            feature = SeqFeature(location, type=feature_type, id=name, strand=1)
            seqrecord.features.append(feature)

    # Add proteins and other features from HXB2
    from operator import attrgetter
    from seqanpy import align_overlap
    from hivwholeseq.genome_info import proteins, chunks
    from hivwholeseq.reference import load_custom_reference
    additional_features_dict = {}
    if 'protein' in additional_features:
        additional_features_dict['protein'] = proteins
    if 'chunk' in additional_features:
        additional_features_dict['chunk'] = chunks

    ref_ann = load_custom_reference('HXB2', 'gb')
    for feagroup, additional_features_grp in additional_features_dict.iteritems():
        for feaname in additional_features_grp:
            if VERBOSE >= 2:
                print feaname,

            fea = ref_ann.features[map(attrgetter('id'), ref_ann.features).index(feaname)]
            seq = fea.extract(ref_ann)
            (score, ali1, ali2) = align_overlap(seqrecord, seq, score_gapopen=-20)
            start = len(ali2) - len(ali2.lstrip('-'))
            end = len(ali2.rstrip('-'))
            end -= ali1[start: end].count('-')

            location = FeatureLocation(start, end)
            if VERBOSE >= 2:
                print 'found:', location

            feature = SeqFeature(location, type=feagroup, id=feaname, strand=1)
            seqrecord.features.append(feature)


def compare_annotations(refseq_new, refseq_old, VERBOSE=0):
    '''Compare annotations of two sequences that are supposed to match'''
    from operator import attrgetter
    feanames_old = map(attrgetter('id'), refseq_old.features)
    feanames_new = map(attrgetter('id'), refseq_new.features)

    feanames = list(feanames_old)
    feanames.extend([fn for fn in feanames_new if fn not in feanames])

    if VERBOSE >= 1:
        print 'Comparison of old and new annotations'
    alert = False
    for fn in feanames:
        if (fn in feanames_old) and (fn not in feanames_new):
            if VERBOSE >= 1:
                print fn, 'present only in old sequence'
            alert = True
            continue
            
        if (fn not in feanames_old) and (fn in feanames_new):
            if VERBOSE >= 1:
                print fn, 'present only in new sequence'
            continue
            
        feaold = refseq_old.features[feanames_old.index(fn)]
        feanew = refseq_new.features[feanames_new.index(fn)]
        
        locold = feaold.location
        locnew = feanew.location

        if str(locold) != str(locnew):
            if VERBOSE >= 1:
                print fn, 'coordinates do not correspond', locold, locnew
            alert = True
            break

    if alert:
        raise ValueError('Not all features are fine.')



# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Annotate initial reference')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--save', action='store_true',
                        help='Save annotated reference to file')
    parser.add_argument('--force', action='store_true',
                        help='Go ahead even if annotations differ from existing sequence')

    args = parser.parse_args()
    VERBOSE = args.verbose
    pnames = args.patients
    use_save = args.save
    use_force = args.force

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[patients.index.isin(pnames)]

    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        if VERBOSE:
            print 'Patient:', patient.name

        fn = patient.get_reference_filename('genomewide')
        refseq = SeqIO.read(fn, 'fasta', alphabet=ambiguous_dna)

        fragment_edges = get_edges_fragments(patient, VERBOSE=VERBOSE)
        annotate_sequence(refseq, VERBOSE=VERBOSE,
                          additional_edges={'fragment': fragment_edges})

        if VERBOSE >= 1:
            for feature in refseq.features:
                if feature.id[0] == 'F':
                    continue

                print feature.id, 
                from hivwholeseq.genome_info import genes
                if feature.type in ('gene', 'protein'):
                    print feature.extract(refseq).seq.translate()
                else:
                    print feature.extract(refseq).seq
            print ''


        try:
            refseq_old = patient.get_reference('genomewide', format='gb')
        except IOError:
            if VERBOSE >= 1:
                print "Old annotated reference not found (that's ok)"
            refseq_old = None

        if refseq_old is not None:
            try:
                compare_annotations(refseq, refseq_old, VERBOSE=VERBOSE)
            except ValueError:
                if use_force:
                    print 'Annotations differ from old sequence'
                else:
                    raise

        if use_save:
            fn_out = patient.get_reference_filename('genomewide', format='gb')
            correct_genbank_features_save(refseq)
            SeqIO.write(refseq, fn_out, 'gb')

