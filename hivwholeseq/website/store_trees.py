# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/11/14
content:    Store phylogenetic trees in a suitable format for the website.
'''
# Modules
import os
import sys
from Bio import Phylo

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.filenames import get_consensi_tree_filename as gfn_in
from hivwholeseq.website.filenames import get_consensi_tree_filename as gfn_out


# Globals
refnames = ['LAI-III']



# Script
if __name__ == '__main__':

    fragments = ['F'+str(i) for i in xrange(1, 7)]
    patients = load_patients()
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        print patient.code, patient.name

        for fragment in fragments:
            fn = patient.get_consensi_tree_filename(fragment)

            if not os.path.isfile(fn):
               continue

            tree = Phylo.read(fn, 'newick')

            # Delete initial reference
            for leaf in tree.get_terminals():
                if 'init' in leaf.name:
                    break
            else:
                raise ValueError('Initial reference not found in tree')
            tree.prune(leaf)

            # Rename nodes without sample names
            for leaf in tree.get_terminals():
                leaf.name = str(int(float(leaf.name.split('_')[0])))+'_days'

            # Ladderize
            tree.ladderize()

            # Write output
            fn_out = gfn_out(patient.code, fragment)
            Phylo.write([tree], fn_out, 'newick')


    print 'All patients'
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
