# vim: fdm=indent
'''
author:     Fabio Zanini
date:       12/08/13
content:    Chop HXB2 into the 6 fragments, which are used as "chromosomes" for
            mapping.
'''
# Modules
import os
import numpy as np
from operator import attrgetter
import Bio.SeqIO as SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

from hivwholeseq.reference import load_custom_reference, save_custom_reference


# Globals
coordinates = {'gene': {'gag': [(789, 2292)],
                        'pol': [(2084, 5096)],
                        'env': [(6224, 8795)],
                        'vif': [(5040, 5619)],
                        'vpr': [(5558, 5850)],
                        'vpu': [(6061, 6310)],
                        'tat': [(5830, 6045), (8378, 8469)],
                        'rev': [(5969, 6045), (8378, 8653)],
                        'nef': [(8796, 9417)],
                       },
               'protein': {'p17': [(789, 1185)],
                           'p24': [(1185, 1878)],
                           'p2': [(1878, 1920)],
                           'p7': [(1920, 2085)],
                           'p1': [(2085, 2133)],
                           'p6': [(2133, 2292)],
                           'PR': [(2252, 2549)],
                           'RT': [(2549, 3869)],
                           'p15': [(3869, 4229)],
                           'IN': [(4229, 5096)],
                           'gp120': [(6314, 7757)],
                           'gp41': [(7757, 8795)],
                          },
               'RNA structure': {"LTR5'": [(0, 634)],
                                 "LTR3'": [(9085, 9719)],
                                 'RRE': None,
                                },
               'other': {'V1': [(6615, 6692)],
                         'V2': [(6693, 6812)],
                         # For V3 Jan wants a specific chunk
                         'V3': None,
                         'V4': [(7377, 7478)],
                         'V5': [(7602, 7634)],
                         'psi': None,
                         'env peptide': None},
               'chunk': {'RT1': [(2549, 2549 + 351)],
                         'RT2': [(2549 + 351,  2549 + 2 * 351)],
                         'RT3': [(2549 + 2 * 351,  2549 + 3 * 351)],
                         'RT4': [(2549 + 3 * 351,  3869)],
                         'IN1': [(4229, 4229 + 351)],
                         'IN2': [(4229 + 351, 4229 + 2 * 351)],
                         'IN3': [(4229 + 2 * 351, 5096)],
                         'gp1201': [(6314, 6314 + 297)],
                        },
              }

def get_coordinates_feature(smat, name):
    '''Get the coordinates of a feature that's missing them'''
    #NOTE: this function is semi-official and does not handle compound features
    from hivwholeseq.genome_info import all_edges, find_region_edges

    edges_chunk = all_edges[name]
    edges = find_region_edges(smat, edges_chunk)
    # Some features must be stripped of primers
    if name in ['V3']:
        edges[0] += len(edges_chunk[0])
        edges[1] -= len(edges_chunk[1])

    return [edges]



# Script
if __name__ == '__main__':

    seqold = load_custom_reference('HXB2', 'gb')

    seqnew = load_custom_reference('HXB2', 'fasta')
    smat = np.array(seqnew)

    print 'Add features'
    for typ, coord_typ in coordinates.iteritems():
        for name, edges in coord_typ.iteritems():
            # If coordinates are missing, grab primers
            if edges is None:
                edges = get_coordinates_feature(smat, name)

            if len(edges) == 1:
                fea = SeqFeature(FeatureLocation(edges[0][0], edges[0][1], strand=+1),
                                 type=typ,
                                 id=name)

            else:
                fea = SeqFeature(CompoundLocation([FeatureLocation(edges[0][0], edges[0][1], strand=+1),
                                                   FeatureLocation(edges[1][0], edges[1][1], strand=+1)]),
                                 type=typ,
                                 id=name)
            
            seqnew.features.append(fea)

    
    print 'Sanity checks'
    for fea in seqnew.features:
        if fea.type in ('gene', 'protein'):
            feaseq = fea.extract(seqnew)

            # vpr has an additional T in HXB2
            if fea.id == 'vpr':
                assert ((len(feaseq) - 1) % 3) == 0
                T_pos = 5771
                prot = (feaseq[: T_pos - fea.location.nofuzzy_start].seq + \
                        feaseq[T_pos + 1 - fea.location.nofuzzy_start:].seq).translate()
            else:
                assert (len(feaseq) % 3) == 0
                prot = feaseq.seq.translate()


            # These genes contain premature stops in HXB2
            if fea.id in ('tat', 'nef'):
                assert prot[-1] == '*'
                assert prot[:-1].count('*') == 1
            elif (fea.type == 'gene') or (fea.id in ('p6', 'IN', 'gp41')):
                assert prot[-1] == '*'
                assert '*' not in prot[:-1]
            else:
                assert '*' not in prot


    print 'Compare to old sequence'
    for fea in seqnew.features:
        if fea.id not in map(attrgetter('id'), seqold.features):
            continue

        feaold = seqold.features[map(attrgetter('id'), seqold.features).index(fea.id)]

        feaseq = fea.extract(seqnew)
        feaoldseq = feaold.extract(seqold)

        # NOTE: add an if check here in case of misannotations
        assert ''.join(feaseq) == ''.join(feaoldseq)


    print 'Importing additional features from old sequence'
    for fea in seqold.features:
        if fea.id not in map(attrgetter('id'), seqnew.features):
            seqnew.features.append(fea)
            print 'Added feature:', fea.id


    print 'Storing new annotated HXB2 sequence'
    save_custom_reference(seqnew, 'HXB2', 'gb')  
