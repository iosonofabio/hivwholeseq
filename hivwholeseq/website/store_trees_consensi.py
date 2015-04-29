# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/11/14
content:    Store phylogenetic trees in a suitable format for the website.

            NOTE: the trees have been JSONed and anonymyzed already.
'''
# Modules
import os
import shutil
import sys
from operator import attrgetter

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.filenames import get_consensi_tree_filename as gfn_in
from hivwholeseq.website.filenames import get_consensi_tree_filename as gfn_out
from hivwholeseq.utils.tree import tree_to_json, correct_minimal_branches
from hivwholeseq.utils.generic import write_json
from hivwholeseq.reference import load_custom_reference



# Globals
refnames = ['LAI-III']



# Script
if __name__ == '__main__':


    # GLOABL CONSENSI TREES
    print 'All patients'
    refseq = load_custom_reference('HXB2', 'gb')
    regions = ([fea.id for fea in refseq.features if len(fea.location.parts) == 1] +
               ['F'+str(i) for i in xrange(1, 7)])

    for region in regions:
        print region,
        fn = gfn_in('all', region, format='json')
        if not os.path.isfile(fn):
            print 'SKIP'
            continue

        # We need to reroot it between the subtypes
        from hivwholeseq.utils.generic import read_json
        from hivwholeseq.utils.tree import tree_from_json, tree_to_json
        tree = tree_from_json(read_json(fn))

        correct_minimal_branches(tree)

        tree.root_at_midpoint()

        tree_json = tree_to_json(tree.root,
                                 fields=('CD4', 'DSI', 'VL',
                                         'confidence',
                                         'muts',
                                         'patient',
                                         'sequence',
                                         'subtype'),
                                )

        # Write output
        fn_out = gfn_out('all', region, format='json')
        write_json(tree_json, fn_out)
        print 'OK'


    # PATIENT CONSENSI TREES
    patients = load_patients()
    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        refseq_gw = patient.get_reference('genomewide', 'gb')
        regions = map(attrgetter('id'), refseq_gw.features) + ['genomewide']

        for region in regions:
            print patient.code, patient.name, region

            fn = patient.get_consensi_tree_filename(region, format='json')
            if not os.path.isfile(fn):
               continue

            tree = patient.get_consensi_tree(region, format='json')

            # Delete initial reference
            for leaf in tree.get_terminals():
                if 'reference' in leaf.name:
                    break
            else:
                raise ValueError('Initial reference not found in tree')
            tree.prune(leaf)

            tree.ladderize()

            correct_minimal_branches(tree)

            # Write output
            fn_out = gfn_out(patient.code, region, format='json')
            tree_json = tree_to_json(tree.root,
                                     fields=('DSI', 'sequence', 'muts',
                                             'VL', 'CD4', 'confidence',
                                             'subtype'),
                                    )
            write_json(tree_json, fn_out, indent=1)