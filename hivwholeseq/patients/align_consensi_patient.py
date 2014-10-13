# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/09/14
content:    Grab all consensi from a patient and align them.
'''
# Modules
import sys
import os
import argparse
import numpy as np
from Bio import SeqIO, AlignIO

from hivwholeseq.sequencing.samples import SampleSeq
from hivwholeseq.patients.patients import load_patients, Patient, SamplePat
from hivwholeseq.tree_utils import build_tree_fasttree


# Functions



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Align consensi',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+', default=['all'],
                        help='Patients to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save alignment to file')
    parser.add_argument('--PCR1', type=int, default=1,
                        help='Take only PCR1 samples [0=off, 1=where both available, 2=always]')
    parser.add_argument('--tree', action='store_true',
                        help='Make phylogenetic tree and save it. Requires --save')
    parser.add_argument('--plot-tree', action='store_true', dest='plot_tree',
                        help='Plot phylogenetic tree. Requires --save and --tree')
    parser.add_argument('--raw', action='store_true',
                        help='Use the raw consensi from premapping')

    args = parser.parse_args()
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    use_save = args.save
    use_PCR1 = args.PCR1
    use_tree = args.tree
    use_plottree = args.plot_tree
    use_raw = args.raw

    patients = load_patients()
    if pnames != ['all']:
        patients = patients.iloc[patients.index.isin(pnames)]

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    alis = {}
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()

        for fragment in fragments:
            if VERBOSE >= 1:
                print pname, fragment,
                if VERBOSE == 1:
                    print ''

            refseq = patient.get_reference(fragment)
            refseq.id = 'reference_'+refseq.id
            refseq.name = 'reference_'+refseq.name
            refseq.description = 'reference '+refseq.description
            seqs = [refseq]

            for i, (samplename, sample) in enumerate(patient.samples.iterrows()):
                if VERBOSE >= 2:
                    print samplename,

                if not use_raw:
                    sample = SamplePat(sample)
                    fn = sample.get_consensus_filename(fragment, PCR=1)
                    if os.path.isfile(fn):
                        cons_seq = SeqIO.read(fn, 'fasta')
                        cons_seq.id = str(patient.times[i])+'_'+cons_seq.id
                        cons_seq.name = cons_seq.id
                        seqs.append(cons_seq)
                        if VERBOSE >= 2:
                            print 'PCR1 OK',

                        if use_PCR1 > 0:
                            if VERBOSE >= 2:
                                print ''
                            continue

                    else:
                        if VERBOSE >= 2:
                            print 'PCR1 not found',

                        if use_PCR1 == 2:
                            if VERBOSE >= 2:
                                print ''
                            continue

                    fn = sample.get_consensus_filename(fragment, PCR=2)
                    if os.path.isfile(fn):
                        cons_seq = SeqIO.read(fn, 'fasta')
                        cons_seq.id = str(patient.times[i])+'_'+cons_seq.id
                        cons_seq.name = cons_seq.id
                        seqs.append(cons_seq)
                        if VERBOSE >= 2:
                            print 'PCR2 OK'
                    else:
                        if VERBOSE >= 2:
                            print 'PCR2 not found'

                else:
                    samples_seq = sample['samples seq']
                    if use_PCR1 == 2:
                        samples_seq = samples_seq.loc[samples_seq.PCR == 1]
                    elif use_PCR1 == 1:
                        if (samples_seq.PCR == 1).any():
                            samples_seq = samples_seq.loc[samples_seq.PCR == 1]

                    if not samples_seq.shape[0]:
                        if VERBOSE >= 2:
                            print 'no sequenced samples found'
                        continue

                    for samplename_seq, sample_seq in samples_seq.iterrows():
                        sample_seq = SampleSeq(sample_seq)
                        fn = sample_seq.get_consensus_filename(fragment)
                        if not os.path.isfile(fn):
                            continue
                        cons_seq = sample_seq.get_consensus(fragment)
                        cons_seq.id = str(patient.times[i])+'_'+cons_seq.id
                        cons_seq.name = cons_seq.id
                        seqs.append(cons_seq)
                        if VERBOSE >= 2:
                            print 'OK'
                        # Check for length and stuff, or collect all
                        break
                    else:
                        if VERBOSE >= 2:
                            print 'no sequenced samples found'
                    
            if VERBOSE >= 2:
                print 'Align',
                sys.stdout.flush()
            from hivwholeseq.mapping_utils import align_muscle
            ali = align_muscle(*seqs, sort=True)
            if VERBOSE >= 2:
                print 'OK'

            if VERBOSE >= 2:
                alim = np.array(ali)
                poly = (alim != alim[0]).any(axis=0)
                print 'Polymorphic sites: ', poly.sum()
                if poly.sum():
                    polyi = poly.nonzero()[0]
                    n_reps = 1 + (poly.sum() - 1) // 50
                    for n_rep in xrange(n_reps):
                        for i in xrange(len(ali)):
                            print '{:>20}'.format(ali[i].id[:20]),
                            print''.join(map(str, alim[i, polyi[n_rep * 50:(n_rep + 1) * 50]]))
                        if n_rep != n_reps - 1:
                            print ''
                
                print ''


            if use_save:
                if VERBOSE >= 2:
                    print 'Save',
                fn_out = patient.get_consensi_alignment_filename(fragment)
                AlignIO.write(ali, fn_out, 'fasta')
                if VERBOSE >= 2:
                    print 'OK'


                if use_tree:
                    if VERBOSE >= 2:
                        print 'Build and save tree',
                        sys.stdout.flush()
                    tree = build_tree_fasttree(fn_out, rootname=seqs[0].id,
                                               VERBOSE=VERBOSE)
                    tree.ladderize()

                    from Bio import Phylo
                    Phylo.write(tree, patient.get_consensi_tree_filename(fragment),
                                'newick')

                    if VERBOSE >= 2:
                        print 'OK'

                    if use_plottree:
                        import matplotlib.pyplot as plt
                        fig, ax = plt.subplots(figsize=(15, 12))
                        Phylo.draw(tree, do_show=False, axes=ax)
                        ax.set_title(pname+', '+fragment)

                        x_max = max(tree.depths().itervalues())
                        ax.set_xlim(0.995, 0.995 + (x_max - 0.995) * 1.4)
                        ax.grid(True)
                        
                        plt.ion()
                        plt.show()





