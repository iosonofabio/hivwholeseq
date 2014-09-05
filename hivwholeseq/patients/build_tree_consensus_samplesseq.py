# vim: fdm=marker
'''
author:     Fabio Zanini
date:       22/06/14
content:    Take all consensi from sequenced samples connected to patient samples,
            make a multiple sequence alignment (MSA), and build a tree with it.

            Seq repetitions and PCR1/2 of the same sample should cluster, then within
            a patient. Else, contamination!
'''
# Modules
import os
import sys
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo

from hivwholeseq.patients.patients import load_patient
from hivwholeseq.patients.patients import load_samples_sequenced as lssp
from hivwholeseq.sequencing.samples import load_samples_sequenced as lss
from hivwholeseq.sequencing.samples import SampleSeq
from hivwholeseq.sequencing.filenames import get_consensus_filename
from hivwholeseq.mapping_utils import align_muscle
from hivwholeseq.patients.filenames import get_multiple_sequence_alignment_sampleseq_filename, \
        get_tree_sampleseq_filename


# Globals


# Functions
def calculate_tree(aln, outgroup=None, VERBOSE=0):
    '''
    takes a biopython alignment and an outgroup sequences, writes a temporary alignment
    and builds a tree using fasttree
    '''
    import subprocess, gzip
    from StringIO import StringIO
    from Bio import SeqIO, Phylo
    
    fasttree_bin = '/ebio/ag-neher/home/fzanini/bin/fasttree'
    home_folder = '/ebio/ag-neher/home/fzanini/'

    tmp_filename = home_folder+'tmp/tmp_'+str(np.random.randint(1000000000))+'.fasta.gz'
    with gzip.open(tmp_filename, 'w') as outfile:
        if outgroup is not None:
            SeqIO.write(outgroup, outfile, 'fasta')
        AlignIO.write(aln, outfile, 'fasta')

    try:
        if VERBOSE >= 2:
            print "building tree of file", tmp_filename
        fasttree_cmd = [fasttree_bin, '-nt']
        ps = subprocess.Popen(('gunzip', '-c', tmp_filename),
                              stdout=subprocess.PIPE)
        fasttree_output = subprocess.check_output(fasttree_cmd,
                                                  stdin=ps.stdout,
                                                  stderr=subprocess.STDOUT)
        ps.wait()

        #all fasttree output is returned as byte string, parse and write tree to file
        biopython_tree = Phylo.read(StringIO(fasttree_output.split('\n')[-2]), 'newick')
        if outgroup is not None:
            # determine the root clade, root, and ladderize
            outgroup_clade = [z for z in biopython_tree.get_terminals()
                              if z.name.startswith(outgroup.name)][0]
            biopython_tree.root_with_outgroup(outgroup_clade)
            biopython_tree.ladderize()
            biopython_tree.prune(outgroup_clade)
            if VERBOSE >= 2:
                print 'Rooted and ladderized'

        biopython_tree.root.branch_length = 0.001

    finally:
        os.remove(tmp_filename)

    return biopython_tree



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Map to initial consensus')
    pats_or_samples = parser.add_mutually_exclusive_group(required=False)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to map (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--fragments', nargs='+',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--no-summary', action='store_false', dest='summary',
                        help='Do not save results in a summary file')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    fragments = args.fragments
    submit = args.submit
    VERBOSE = args.verbose
    summary = args.summary


    samples_pat = lssp()

    # Collect all sequenced samples from patients
    if pnames is not None:
        samples_seq = []
        for pname in pnames:
            patient = load_patient(pname)
            patient.discard_nonsequenced_samples()
            for samplename_pat, sample_pat in patient.samples.iterrows():
                samples_seq.append(sample_pat['samples seq'])
        samples_seq = pd.concat(samples_seq)

    else:
        samples_seq = lss()
        if samplenames is not None:
            ind = samples_pat.index.isin(samplenames)
            if ind.sum():
                samplenames_pat = samples_pat.index[ind]
                samples_seq = samples_seq.loc[samples_seq['patient sample'].isin(samplenames_pat)]
            else:
                samples_seq = samples_seq.loc[samples_seq.index.isin(samplenames)]
        else:
            samples_seq = samples_seq.loc[samples_seq['patient sample'] != 'nan']

    if VERBOSE >= 2:
        print 'samples', samples_seq.index.tolist()

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    for fragment in fragments:
        if VERBOSE >= 1:
            print fragment

        consensi = []
        for samplename, sample in samples_seq.iterrows():
            sample = SampleSeq(sample)

            samplename_pat = sample['patient sample']
            sample_pat = samples_pat.loc[samplename_pat] 
            sample['patient'] = pname = sample_pat.patient
            date = sample_pat.date
            PCR = int(sample.PCR)
            data_folder = sample.sequencing_run.folder
            seq_run = sample['seq run']
            adaID = sample.adapter

            if VERBOSE >= 1:
                print samplename, seq_run, adaID,

            cons_fn = get_consensus_filename(data_folder, adaID, fragment)
            if os.path.isfile(cons_fn):
                if VERBOSE >= 3:
                    print 'consensus found'
                elif VERBOSE >= 1:
                    print ''
                consensus = SeqIO.read(cons_fn, 'fasta')
                consensus.description = consensus.description+', patient '+sample.patient
                consensus.id = consensus.name = pname+'-'+str(date.date())+'-'+samplename
                consensi.append(consensus)
            else:
                print 'ERROR: consensus file not found'

        if VERBOSE:
            print 'N consensi:', len(consensi)
                
        # Align
        if VERBOSE:
            print 'Align'
        ali = align_muscle(*consensi)
        ali_fn = get_multiple_sequence_alignment_sampleseq_filename(fragment)
        AlignIO.write(ali, ali_fn, 'fasta')

        # Make tree
        if VERBOSE:
            print 'Build tree'
        tree_fn = get_tree_sampleseq_filename(fragment)
        tree = calculate_tree(ali, VERBOSE=VERBOSE)
        Phylo.write(tree, tree_fn, 'newick')
