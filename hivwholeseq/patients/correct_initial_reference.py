# vim: fdm=marker
'''
author:     Fabio Zanini
date:       13/11/14
content:    Correct initial reference to look more like the initial consensus,
            without changing the coordinates, using allele frequencies.
'''
# Modules
import os
import sys
import stat
import argparse
import shutil
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio import SeqIO

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.patients import load_patients, Patient



# Functions
def save_protect(fn, fn_old, VERBOSE=0):
    '''Save new refernce and protect old one'''
    # Old files get protected, check for it
    if os.path.isfile(fn_old):
        is_writable = bool(os.stat(fn_old).st_mode & stat.S_IWUSR)
    else:
        is_writable = True

    if not is_writable:
        print 'Old reference is protected, cannot save'
    else:
        shutil.copy(fn, fn_old)
        SeqIO.write(seq, fn, 'fasta')
        # Protect old reference (this can be done once only!)
        os.chmod(fn_old, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
        if VERBOSE >= 1:
            print 'New reference saved and protected'


def merge_sequences_fragments(seqs, VERBOSE=0):
    '''Merge sequences from consecutive fragments'''
    from seqanpy import align_ladder

    seqs = map(''.join, seqs)
    seq = [seqs[0]]
    for seq2 in seqs[1:]:
        seq1 = seq[-1]
        (score, ali1, ali2) = align_ladder(seq1, seq2, score_gapopen=-20)
        start2 = len(ali2) - len(ali2.lstrip('-'))
        end1 = len(ali1.rstrip('-'))
        len_overlap = end1 - start2

        # Trust the first sequence in the first half, the second in the second
        overlap = ali1[start2: start2 + len_overlap / 2] + \
                  ali2[start2 + len_overlap / 2: end1]

        seq[-1] = ali1[:start2]
        seq.append(overlap)
        seq.append(ali2[end1:])

    seq = ''.join(seq)
    return seq



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get time to boundary (fix/loss)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 genomewide)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Write new reference')
    parser.add_argument('--recover', action='store_true',
                        help='Recover old reference')

    args = parser.parse_args()
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    use_save = args.save
    use_recover = args.recover

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    if use_recover:
        for pname, patient in patients.iterrows():
            print pname
            patient = Patient(patient)
            patient.discard_nonsequenced_samples()

            for fragment in fragments:
                if VERBOSE >= 1:
                    print fragment

                fn = patient.get_reference_filename(fragment)
                fn_old = fn.replace('.fasta', '_old.fasta')
                if not os.path.isfile(fn_old):
                    print 'Old reference not found, skipping'
                    continue
                shutil.copy(fn_old, fn)
                os.chmod(fn, (stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH | \
                              stat.S_IWUSR))
                os.chmod(fn_old, (stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH | \
                                  stat.S_IWUSR))
                os.remove(fn_old)
                if VERBOSE >= 1:
                    print 'Old reference recovered and protection removed'
                continue
        sys.exit()

    if 'genomewide' in fragments:
        do_genomewide = True
        del fragments[fragments.index('genomewide')]
    else:
        do_genomewide = False

    for pname, patient in patients.iterrows():
        print pname
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()

        for fragment in fragments:
            if VERBOSE >= 1:
                print fragment

            ref = patient.get_reference(fragment)
            refm = np.array(ref, 'S1')
            af0 = patient.get_initial_allele_frequencies(fragment)
            consm = alpha[af0.argmax(axis=0)]

            seqm = consm.copy()
            # Gaps in a reference are not wished for
            if '-' in seqm:
                seqm[seqm == '-'] = refm[seqm == '-']

            seq = SeqRecord(Seq(''.join(seqm), ambiguous_dna),
                            id=ref.id,
                            name=ref.name,
                            description=ref.description)

            if VERBOSE >= 2:
                ind = (refm != seqm).nonzero()[0]
                if not len(ind):
                    print 'Nothing changed'
                else:
                    'Position changed:',
                    for i in ind:
                        print i, refm[i], '->', seqm[i]

            if use_save:
                fn = patient.get_reference_filename(fragment, 'fasta')
                fn_old = fn.replace('.fasta', '_old.fasta')
                save_protect(fn, fn_old, VERBOBE=VERBOSE)


        if do_genomewide:
            seqs = [patient.get_reference('F'+str(i)) for i in xrange(1, 7)]
            seq = merge_sequences_fragments(seqs, VERBOSE=VERBOSE)

            seq = SeqRecord(Seq(seq, ambiguous_dna),
                            id=pname+'_genomewide',
                            name=pname+'_genomewide',
                            description='Genomewide reference for patient '+pname)

            ref = patient.get_reference('genomewide')

            if VERBOSE >= 2:
                from seqanpy import align_global
                from hivwholeseq.sequence_utils import pretty_print_pairwise_ali
                (score, ali1, ali2) = align_global(ref, seq, score_gapopen=-20)
                pretty_print_pairwise_ali((ali1, ali2),
                                          name1='Old ref', name2='New ref',
                                          width=100)

            # TODO: resplit sequences to make sure we cover the whole F5a, F3c,
            # etc. THIS CHANGES THE COORDINATES!

            if use_save:
                fn = patient.get_reference_filename('genomewide', 'fasta')
                fn_old = fn.replace('.fasta', '_old.fasta')
                save_protect(fn, fn_old, VERBOSE=VERBOSE)


