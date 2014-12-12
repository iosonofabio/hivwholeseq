#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       03/12/14
content:    Compute alignments of haplotypes in a few regions and store them for
            the website.
'''
# Modules
import os
import sys
import argparse
import numpy as np
import StringIO
from Bio import AlignIO
from Bio import Phylo

from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.sequence_utils import build_msa_haplotypes as build_msa
from hivwholeseq.website.filenames import get_precompiled_alignments_filename
from hivwholeseq.fork_cluster import fork_store_haplotypes_website as fork_self
from hivwholeseq.sequence_utils import align_muscle
from hivwholeseq.tree_utils import build_tree_fasttree
from hivwholeseq.patients.get_local_trees import get_region_count_trajectories
from hivwholeseq.patients.filenames import get_consensi_tree_filename as gfn_in
from hivwholeseq.website.filenames import get_consensi_tree_filename as gfn_out
from hivwholeseq.website.filenames import get_consensi_alignment_filename


# Globals
regions_all = ['PR', 'V3', 'psi']
freqmin = 0.01



# Functions
def store_alignments(alis, filename):
    '''Store alignments to file'''
    import zipfile, zlib
    with zipfile.ZipFile(filename, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
        for i, ali in enumerate(alis):
            f = StringIO.StringIO()
            AlignIO.write(ali['ali'], f, 'fasta')
            zf.writestr(str(ali['time'])+'_days.fasta', f.getvalue())


def annotate_tree_consensi(tree, times):
    '''Change leaf names to include more info'''
    # Reroot
    iseqroot = 0
    for leaf in tree.get_terminals():
        iseq = int(leaf.name[3:]) - 1
        if iseq == iseqroot:
            break
    else:
        raise ValueError('Not able to reroot!')
    tree.root_with_outgroup(leaf)

    # Change name into the time
    for leaf in tree.get_terminals(): 
        iseq = int(leaf.name[3:]) - 1
        leaf.name = str(times[iseq])+'_days'


def annotate_tree_minor(tree, hft, times, seqs):
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

    # Parse input args
    parser = argparse.ArgumentParser(description='Perform PCA on the data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--regions', nargs='+', default=regions_all,
                        help='Regions to store (e.g. V3 PR)')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    submit = args.submit

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    if submit:
        for pname, patient in patients.iterrows():
            for region in regions:
                fork_self(pname, region, VERBOSE=1)
        sys.exit()

    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        for region in regions:
            print patient.code, patient.name, region

            # Compute alignments
            print 'Get haplotypes'
            haplos = patient.get_local_haplotype_trajectories(region, 0, '+oo')
            print 'Align'
            alis = [{'time': patient.times[it],
                     'ali': build_msa(h)}
                    for h, it in zip(*haplos)]

            # Write to file
            print 'Write output file'
            fn_out = get_precompiled_alignments_filename(patient.code, region)
            mkdirs(os.path.dirname(fn_out))
            store_alignments(alis, fn_out)

            # Trees (we make them de novo from extant alignments)
            print patient.code, patient.name, region

            print 'Get haplotypes'
            # NOTE: this function is kinda ridiculous, since we have "alis", but
            # it makes the code more modular so it's staying for now
            (hct, ind, seqs) = get_region_count_trajectories(patient, region,
                                                             VERBOSE=2)
            hct = hct.T
            hft = 1.0 * hct / hct.sum(axis=0)
        
            # Exclude too rare haplos
            indseq = (hft >= freqmin).any(axis=1)
            seqs = seqs[indseq]
            hft = hft[indseq]
        
            times = patient.times[ind]
        
            # 1. Make consensi tree
            consensi = seqs[hft.argmax(axis=0)]
            
            print 'Align consensi'
            consali = align_muscle(*consensi, sort=True)

            print 'Build local tree of consensi'
            tree = build_tree_fasttree(consali, VERBOSE=2)

            # Change names and so on
            annotate_tree_consensi(tree, times)

            # Ladderize in place
            tree.ladderize()

            # Write output
            fn_out = gfn_out(patient.code, region)
            Phylo.write([tree], fn_out, 'newick')

            # 2. Make minor haplotype tree
            print 'Align sequences'
            seqsali = align_muscle(*seqs, sort=True)
        
            print 'Build local tree'
            tree = build_tree_fasttree(seqsali, VERBOSE=2)

            # Change names and so on
            annotate_tree_minor(tree, hft, times, seqs)

            # Ladderize in place
            tree.ladderize()

            # Write output
            fn_out = gfn_out(patient.code, region+'minor')
            Phylo.write([tree], fn_out, 'newick')

            # 3. Save consensi
            print 'Save consensi'
            for i, seq in enumerate(consali):
                tstr = str(times[i])
                consali[i].id = '_'.join([patient.code, region, tstr, 'days'])
                consali[i].name = consali[i].id
                consali[i].description = ', '.join([patient.code, region, tstr+' days'])
            
            # Write output
            fn_out = get_consensi_alignment_filename(patient.code, region,
                                                     format='fasta')
            AlignIO.write(consali, fn_out, 'fasta')


