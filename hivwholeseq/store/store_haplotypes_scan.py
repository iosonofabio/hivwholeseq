#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/12/14
content:    Store local haplotypes.
'''
# Modules
import os
import sys
import argparse
from operator import itemgetter, attrgetter
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from Bio import Phylo

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.utils.mapping import align_muscle
from hivwholeseq.utils.argparse import PatientsAction
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.utils.exceptions import RoiError
from hivwholeseq.utils.tree import build_tree_fasttree
from hivwholeseq.store.store_tree_consensi import annotate_tree
from hivwholeseq.utils.nehercook.ancestral import ancestral_sequences
from hivwholeseq.utils.tree import tree_to_json
from hivwholeseq.utils.generic import write_json
from hivwholeseq.store.store_haplotypes_alignment_tree import (
    expand_annotate_alignment,
    annotate_tree_time_freq_count, plot_tree)
from hivwholeseq.cluster.fork_cluster import fork_store_haplotypes_scan as fork_self



# Functions



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Store local haplotypes',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', action=PatientsAction,
                        help='Patients to analyze')
    parser.add_argument('--width', type=int, default=400,
                        help='Width of the sliding window [bp]')
    parser.add_argument('--gap', type=int, default=100,
                        help='Gap between two windows [bp]')
    parser.add_argument('--start', type=int, default=0,
                        help='Start position in genome [bp]')
    parser.add_argument('--end', type=int, default=10000,
                        help='End position in genome [bp]')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--freqmin', type=float, default=0.01,
                        help='Minimal frequency to keep the haplotype')
    parser.add_argument('--countmin', type=int, default=3,
                        help='Minimal observations to keep the haplotype')
    parser.add_argument('--plot', action='store_true',
                        help='Plot local haplotype trajectories')
    parser.add_argument('--submit', action='store_true',
                        help='Submit to the cluster')
    parser.add_argument('--save', action='store_true',
                        help='Save alignment to file')

    args = parser.parse_args()
    pnames = args.patients
    width = args.width
    gap = args.gap
    start = args.start
    end = args.end
    VERBOSE = args.verbose
    freqmin = args.freqmin
    countmin = args.countmin
    submit = args.submit
    use_plot = args.plot
    use_save = args.save

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    data = []
    for pname, patient in patients.iterrows():
        if VERBOSE >= 1:
            print patient.code, start, end

        if submit:
            fork_self(patient.code, width, gap, start, end, VERBOSE=VERBOSE,
                      freqmin=freqmin, countmin=countmin)
            continue

        patient = Patient(patient)
        ref = patient.get_reference('genomewide')
        L = len(ref)

        win_start = start
        while win_start + width - gap < min(L, end):
            win_end = min(win_start + width, end, L)

            if VERBOSE >= 1:
                print patient.code, win_start, win_end
    
            if VERBOSE >= 2:
                print 'Get region haplotypes'
            try:
                datum = patient.get_local_haplotype_count_trajectories(\
                               'genomewide',
                               start=win_start,
                               end=win_end,
                               filters=['noN',
                                        'mincount='+str(countmin),
                                        'freqmin='+str(freqmin),
                                       ],
                               VERBOSE=VERBOSE,
                               align=True,
                               return_dict=True)
            except RoiError:
                win_start += gap
                continue

            if not len(datum['ind']):
                win_start += gap
                continue

            datum['times'] = patient.times[datum['ind']]
            datum['pcode'] = patient.code
            datum['window'] = (win_start, win_end)
            data.append(datum)

            if use_save:
                if VERBOSE >= 2:
                    print 'Save to file'
                
                rname = 'scan_'+str(win_start)+'-'+str(win_end)
                fn_out = patient.get_haplotype_count_trajectory_filename(rname)
                mkdirs(os.path.dirname(fn_out))
                np.savez_compressed(fn_out,
                                    hct=datum['hct'],
                                    ind=datum['ind'],
                                    times=datum['times'],
                                    seqs=datum['seqs'],
                                    ali=datum['alim'],
                                   )

            if VERBOSE >= 2:
                print 'Build tree'
            times = datum['times']
            alim = datum['alim']
            hct = datum['hct']
            hft = 1.0 * hct / hct.sum(axis=0)
            ali = expand_annotate_alignment(alim, hft, hct, times,
                                            freqmin=freqmin,
                                            VERBOSE=VERBOSE)
            tree = build_tree_fasttree(ali, VERBOSE=VERBOSE)

            if VERBOSE >= 2:
                print 'Infer ancestral sequences'
            a = ancestral_sequences(tree, ali, alphabet='ACGT-N', copy_tree=False,
                                    attrname='sequence', seqtype='str')
            a.calc_ancestral_sequences()
            a.cleanup_tree()

            if VERBOSE >= 2:
                print 'Annotate tree'
            annotate_tree_time_freq_count(tree, ali)
            annotate_tree(patient, tree, VERBOSE=VERBOSE)

            if VERBOSE >= 2:
                print 'Ladderize tree'
            tree.ladderize()

            if use_save:
                if VERBOSE >= 2:
                    print 'Save tree (JSON)'
                rname = 'scan_'+str(win_start)+'-'+str(win_end)
                fn = patient.get_local_tree_filename(rname, format='json')
                mkdirs(os.path.dirname(fn))
                tree_json = tree_to_json(tree.root,
                                         fields=('DSI', 'sequence',
                                                 'muts',
                                                 'VL', 'CD4',
                                                 'frequency',
                                                 'count',
                                                 'confidence'),
                                        )
                write_json(tree_json, fn)

            if use_plot:
                if VERBOSE >= 2:
                    print 'Plot'
                plot_tree(tree, title=patient.code+', '+str(win_start)+'-'+str(win_end))

            win_start += gap
