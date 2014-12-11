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

from hivwholeseq.sequence_utils import align_muscle
from hivwholeseq.tree_utils import build_tree_fasttree
from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.get_local_trees import get_region_count_trajectories
from hivwholeseq.patients.filenames import get_consensi_tree_filename as gfn_in
from hivwholeseq.website.filenames import get_consensi_tree_filename as gfn_out


# Globals
refnames = ['LAI-III']
regions = ['V3', 'psi', 'PR']
freqmin = 0.01


# Functions
def annotate_tree(tree, hft, times, seqs):
    '''Change leaf names to include more info'''
    # Reroot
    iseqroot = hft[:, 0].argmax()
    for leaf in tree.get_terminals():
        iseq = int(leaf.name[3:]) - 1
        if iseq == iseqroot:
            break
    else:
        raise ValueError('Not able to reroot!')
    tree.root_with_outgroup(leaf)

    # Annotate leaves
    for leaf in tree.get_terminals():
        iseq = int(leaf.name[3:]) - 1
        it = hft[iseq].argmax()
        leaf.day = times[it]
        leaf.hf = hft[iseq, it]
        leaf.seq = seqs[iseq]
        leaf.name = (leaf.seq+'_'+str(int(leaf.day))+'days_fmax_'+
                     '{:1.2f}'.format(leaf.hf))



# Script
if __name__ == '__main__':

    fragments = ['F'+str(i) for i in xrange(1, 7)]
    patients = load_patients()

    # Patient by patient
    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        # Fragment trees (we recycle extant files for speed reasons)
        for fragment in fragments:
            print patient.code, patient.name, fragment

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

            # Ladderize in place
            tree.ladderize()

            # Write output
            fn_out = gfn_out(patient.code, fragment)
            Phylo.write([tree], fn_out, 'newick')


        # Region trees (we make them de novo from extant alignments)
        for region in regions:
            print patient.code, patient.name, region

            print 'Get haplotypes'
            (hct, ind, seqs) = get_region_count_trajectories(patient, region,
                                                             VERBOSE=2)
            hct = hct.T
            hft = 1.0 * hct / hct.sum(axis=0)
        
            # Exclude too rare haplos
            indseq = (hft >= freqmin).any(axis=1)
            seqs = seqs[indseq]
            hft = hft[indseq]
        
            times = patient.times[ind]
        
            print 'Align sequences'
            seqsali = align_muscle(*seqs, sort=True)
        
            print 'Build local tree'
            tree = build_tree_fasttree(seqsali, VERBOSE=2)

            # Change names and so on
            annotate_tree(tree, hft, times, seqs)

            # Ladderize in place
            tree.ladderize()

            # Write output
            fn_out = gfn_out(patient.code, region)
            Phylo.write([tree], fn_out, 'newick')


    # Global trees
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
