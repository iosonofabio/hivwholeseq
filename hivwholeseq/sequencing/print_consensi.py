# vim: fdm=marker
'''
author:     Fabio Zanini
date:       30/06/14
content:    Print a table of all consensi with some summary description.
'''
# Modules
import os
import numpy as np
import pandas as pd
import argparse
from Bio import SeqIO

from hivwholeseq.sequencing.samples import SampleSeq
from hivwholeseq.patients.patients import load_patients, load_patient, Patient
from hivwholeseq.patients.patients import load_samples_sequenced as lssp
from hivwholeseq.sequencing.samples import load_samples_sequenced as lss
from hivwholeseq.sequencing.filenames import get_consensus_filename, get_reference_premap_filename
from seqanpy import align_global


# Functions
def score_consensus(sample, VERBOSE=0):
    '''Score a consensus based on completeness and quality'''
    data_folder = sample.sequencing_run.folder
    adaID = sample.adapter

    frag_spec = filter(lambda x: fragment in x, sample.regions_complete)
    if not len(frag_spec):
        field = ''
        return (True, '')

    fn = get_consensus_filename(data_folder, adaID, fragment)
    if not os.path.isfile(fn):
        return (False, 'MISS')

    frag_spec = frag_spec[0]
    fn_ref = get_reference_premap_filename(data_folder, adaID, frag_spec)
    if not os.path.isfile(fn_ref):
        if frag_spec[:3] == 'F3a':
            frag_spec = frag_spec.replace('a', '')
            fn_ref = get_reference_premap_filename(data_folder, adaID, frag_spec)
            if not os.path.isfile(fn_ref):
                return (False, 'MISSREF')
        else:
            return (False, 'MISSREF')

    ref = SeqIO.read(fn_ref, 'fasta')
    cons = SeqIO.read(fn, 'fasta')
    if len(cons) < len(ref) - 200:
        return (False, 'SHORT')
    elif len(cons) > len(ref) + 200:
        return (False, 'LONG')

    #ali = align_global(str(ref.seq), str(cons.seq), band=200)
    #alim1 = np.fromstring(ali[1], 'S1')
    #alim2 = np.fromstring(ali[2], 'S1')
    #if (alim1 != alim2).sum() >
    return (True, 'OK')



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Map to initial consensus')
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument('--runs', nargs='+',
                        help='Sequencing runs to analyze (e.g. Tue28)')
    target.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    target.add_argument('--samples', nargs='+',
                        help='Samples to map')
    parser.add_argument('--fragments', nargs='+',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of read pairs to map (for testing)')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads to use for mapping')
    parser.add_argument('--skiphash', action='store_true',
                        help=argparse.SUPPRESS)
    parser.add_argument('--no-summary', action='store_false', dest='summary',
                        help='Do not save results in a summary file')
    parser.add_argument('--chunks', type=int, nargs='+', default=[None],
                        help='Only map some chunks (cluster optimization): 0 for automatic detection')

    args = parser.parse_args()
    seq_runs = args.runs
    pnames = args.patients
    samplenames = args.samples
    fragments = args.fragments
    submit = args.submit
    VERBOSE = args.verbose
    threads = args.threads
    n_pairs = args.maxreads
    skip_hash = args.skiphash
    summary = args.summary
    only_chunks = args.chunks

    if seq_runs is not None:
        samples_seq = lss()
        samples_seq = samples_seq.loc[samples_seq['seq run'].isin(seq_runs)]

    elif pnames is not None:
        samples_pat = lssp()
        samples_seq = []
        for pname in pnames:
            patient = load_patient(pname)
            patient.discard_nonsequenced_samples()
            for samplename_pat, sample_pat in patient.samples.iterrows():
                samples_seq.append(sample_pat['samples seq'])
        samples_seq = pd.concat(samples_seq)

    else:
        samples_pat = lssp()
        samples_seq = lss()
        ind = samples_pat.index.isin(samplenames)
        if ind.sum():
            samplenames_pat = samples_pat.index[ind]
            samples_seq = samples_seq.loc[samples_seq['patient sample'].isin(samplenames_pat)]
        else:
            samples_seq = samples_seq.loc[samples_seq.index.isin(samplenames)]

    if VERBOSE >= 2:
        print 'samples', samples_seq.index.tolist()

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    linelen = 49 + 11 * len(fragments)
    print '-' * linelen
    for samplename, sample in samples_seq.iterrows():
        sample = SampleSeq(sample)
        seq_run = sample['seq run']
        data_folder = sample.sequencing_run.folder
        adaID = sample.adapter

        line = ' | '.join(['{:<27s}'.format(sample.name),
                           '{:<8s}'.format(seq_run),
                           '{:<6s}'.format(adaID),
                           '',
                          ])

        for fragment in fragments:
            (score, field) = score_consensus(sample, VERBOSE=VERBOSE)

            # NOTE: Tue52 is essentially failed, maybe we discard
            if (seq_run == 'Tue52') and len(field):
                field = '('+field+')'

            line = line+'{:^8s}'.format(field)+' | '

        print line[:-1]

    print '-' * linelen
            


