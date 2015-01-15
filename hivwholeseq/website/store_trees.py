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

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.filenames import get_consensi_tree_filename as gfn_in
from hivwholeseq.website.filenames import get_consensi_tree_filename as gfn_out
from hivwholeseq.tree_utils import tree_to_json
from hivwholeseq.generic_utils import write_json



# Globals
refnames = ['LAI-III']



# Script
if __name__ == '__main__':

    patients = load_patients()

    # Patient by patient
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

            # Write output
            fn_out = gfn_out(patient.code, region, format='json')
            tree_json = tree_to_json(tree.root,
                                     fields=('DSI', 'sequence', 'muts',
                                             'VL', 'CD4', 'confidence'),
                                    )
            write_json(tree_json, fn_out, indent=1)

    # FIXME: make json global trees
    import sys; sys.exit()

    # Global trees
    print 'All patients'
    # TODO: make alignments of other regions
    fragments = ['F'+str(i) for i in xrange(1, 7)]
    for fragment in fragments:
        print fragment
        fn = gfn_in('all', fragment)
        if not os.path.isfile(fn):
           continue

        tree = Phylo.read(fn, 'newick')

        for pname, patient in patients.iterrows():
            patient = Patient(patient)
            print patient.code, patient.name

            for leaf in tree.get_terminals():
                if patient.name not in leaf.name:
                    continue

                found = False
                leaf_sn = leaf.name.split('_')[-1]
                for samplename, sample in patient.samples.iterrows():
                    if samplename == leaf_sn:
                        found = True
                        break

                if found:
                    label = patient.code+'_'+str(sample['days since infection'])+'_days'
                    leaf.name = label
        
        # Prune wrong samples and unknown references
        leaves_to_prune = []
        for leaf in tree.get_terminals():
            if leaf.name[0] == 'p':
                continue

            for refname in refnames:
                if refname in leaf.name:
                    leaf.name = refname
                else:
                    leaves_to_prune.append(leaf)

        for leaf in leaves_to_prune:
            tree.prune(leaf)

        # Write output
        fn_out = gfn_out('all', fragment)
        Phylo.write([tree], fn_out, 'newick')
