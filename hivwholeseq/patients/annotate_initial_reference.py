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
from hivwholeseq.sequence_utils import correct_genbank_features_save




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


def annotate_sequence(seqrecord, additional_edges={}, VERBOSE=0):
    '''Annotate a consensus with the genes and stuff (in place)'''
    # TODO: what do we do with genes that do not start/end where they are
    # supposed to? Do we follow biology and track their new locations?
    from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
    from hivwholeseq.genome_info import gene_edges, RNA_structure_edges, \
            other_edges, find_region_edges, find_region_edges_multiple, \
            locate_gene
    edge_dict = {'genes': gene_edges,
                 'RNA structures': RNA_structure_edges,
                 'other': other_edges}
    edge_dict.update(additional_edges)

    if VERBOSE:
        print 'Features:', ', '.join(edge_dict.iterkeys())

    smat = np.array(seqrecord)

    for feature_type in edge_dict:
        edges_all = edge_dict[feature_type]
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

    args = parser.parse_args()
    VERBOSE = args.verbose
    pnames = args.patients
    use_save = args.save

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
                          additional_edges={'fragments': fragment_edges})

        if VERBOSE >= 1:
            for feature in refseq.features:
                if feature.id[0] == 'F':
                    continue

                print feature.id, 
                from hivwholeseq.genome_info import genes
                if feature.id in genes:
                    print feature.extract(refseq).seq.translate()
                else:
                    print feature.extract(refseq).seq
            print ''


        if use_save:
            fn_out = patient.get_reference_filename('genomewide', format='gb')
            correct_genbank_features_save(refseq)
            SeqIO.write(refseq, fn_out, 'gb')

        # TODO: propagate to the single fragment references now
