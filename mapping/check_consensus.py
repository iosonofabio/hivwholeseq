# vim: fdm=marker
'''
author:     Fabio Zanini
date:       24/09/13
content:    Check iterative consensus with reference (if appropriate) and between
            fragments at overlaps.
'''
# Modules
import argparse
import numpy as np
import Bio.SeqIO as SeqIO
import Bio.AlignIO as AlignIO

from mapping.datasets import MiSeq_runs
from mapping.primer_info import primers_inner as pri
from mapping.adapter_info import load_adapter_table
from mapping.filenames import get_consensus_filename



# Globals
# Reference samples for various MiSeq runs
references = {28: {2: 'NL4-3', 4: 'SHIV_SF162'}}



# Functions
def find_fragment(refseq, seq):
    '''Find a fragment location in a sequence (do not align)'''
    refs = str(refseq.seq)
    s = str(seq.seq)

    # Align the seq back to the reference
    for bsize in [40, 30, 20, 15, 10, 7]:
        i = 0
        while i + bsize < len(s):
            start = refs.find(s[i:i+bsize])
            if start != -1:
                start -= i
                break
            i += bsize
        if start != -1:
            break
    if start == -1:
        return None

    for bsize in [40, 30, 20, 15, 10, 7]:
        i = len(seq)
        while i - bsize > 0:
            end = refs[start:].find(s[i - bsize: i])
            if end != -1:
                end += bsize + len(seq) - i + start
                break
            i -= bsize
        if end != -1:
            break

    refsf = refseq[start:end]
    return refsf


def find_overlap(seq1, seq2):
    '''Find the overlap between two consecutive fragments'''
    from Bio import pairwise2

    s1 = str(seq1.seq)
    s2 = str(seq2.seq)

    ali1 = pairwise2.align.localms(s1[-400:], s2[:15], 2, -1, -2.5, -0.5)
    start = len(s1) - len(ali1[0][1].lstrip('-'))

    ali2 = pairwise2.align.localms(s1[-15:], s2[:400], 2, -1, -2.5, -0.5)
    end = len(ali2[0][0].rstrip('-'))

    return (seq1[start:], seq2[:end])


def align_muscle(refseq, seq):
    '''Align two sequences via MUSCLE'''
    from Bio.Align.Applications import MuscleCommandline
    muscle_cline = MuscleCommandline(diags=True)
    import subprocess as sp
    child = sp.Popen(str(muscle_cline),
                     stdin=sp.PIPE,
                     stdout=sp.PIPE,
                     stderr=sp.PIPE,
                     shell=True)
    SeqIO.write([refseq, seq], child.stdin, "fasta")
    child.stdin.close()
    child.stderr.close()
    align = AlignIO.read(child.stdout, "fasta")
    child.stdout.close()
    return align



# Script
if __name__ == '__main__':

    # Parse input args: this is used to call itself recursively
    parser = argparse.ArgumentParser(description='Map HIV reads recursively')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--reference', default='',
                        help='Also compare to a reference (e.g. NL4-3)')

    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    reference = args.reference

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = load_adapter_table(data_folder)['ID']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    for adaID in adaIDs:

        # If there is a reference, check it
        if (miseq_run in references.keys()) and \
           (adaID in references[miseq_run]) and (not reference):
            reference = references[miseq_run][adaID]
        if reference:
            refseq = SeqIO.read(data_folder+'reference/'+reference+'.fasta', 'fasta')
            
            for fragment in fragments:
                seq = SeqIO.read(get_consensus_filename(data_folder, adaID, fragment),
                                 'fasta')
                refsf = find_fragment(refseq, seq)

                # Align
                align = align_via_muscle(refsf, seq)

                # Look for mutations
                ali = np.array(align)
                muts = (ali[1] != ali[0]).nonzero()[0]
                if len(muts):
                    print adaID, fragment, muts

                #FIXME
                #import ipdb; ipdb.set_trace()

        # Check the number of Ns (did we cover enough with 1000 reads?)
        for fragment in fragments:
            seq = SeqIO.read(get_consensus_filename(data_folder, adaID, fragment),
                             'fasta')
            seqmat = np.array(seq)
            n_N = (seqmat == 'N').sum()
            if (n_N > 30) or (VERBOSE >= 3):
                print 'Number of Ns:', n_N

        # Check the primers
        for fragment in fragments:
            if fragment == 'F5':
                frag_prim = dataset['primerF5'][dataset['adapters'].index(adaID)]
            else:
                frag_prim = fragment
    
            (prim_fwd, prim_rev) = pri[frag_prim]

            from Bio import pairwise2
            seq = SeqIO.read(get_consensus_filename(data_folder, adaID, fragment),
                             'fasta')

            ali = pairwise2.align.localms(prim_fwd, str(seq.seq)[:len(prim_fwd)],
                                          2, -1, -2.5, -0.5)
            if VERBOSE >= 3:
                print 'fwd primer'
                print ali
            alinp = np.array(map(list, ali[0][:2]))
            alinp = alinp[:, -len(ali[0][1].strip('-')):len(ali[0][0].rstrip('-'))]
            frac_diff = (alinp[0] != alinp[1]).mean()
            if (frac_diff > 0.2) or (VERBOSE >= 3):
                print 'Fwd primer:', adaID, fragment, frac_diff

            ali = pairwise2.align.localms(prim_rev, str(seq.seq)[-len(prim_rev):],
                                          2, -1, -2.5, -0.5)
            if VERBOSE >= 3:
                print 'rev primer'
                print ali
            alinp = np.array(map(list, ali[0][:2]))
            alinp = alinp[:, -len(ali[0][0].strip('-')):len(ali[0][1].rstrip('-'))]
            frac_diff = (alinp[0] != alinp[1]).mean()
            if (frac_diff > 0.2) or (VERBOSE >= 3):
                print 'Rev primer:', adaID, fragment, frac_diff

        # Look at the overlapping regions
        for fragment1 in fragments:
            fragment2 = 'F'+str(int(fragment1[-1])+1)
            if fragment2 not in fragments:
                continue

            seq1 = SeqIO.read(get_consensus_filename(data_folder, adaID, fragment1),
                              'fasta')
            seq2 = SeqIO.read(get_consensus_filename(data_folder, adaID, fragment2),
                              'fasta')

            # Find overlap
            (seqo1, seqo2) = find_overlap(seq1, seq2)

            # Align
            align = align_muscle(seqo1, seqo2)

            # Look for mutations
            ali = np.array(align)
            muts = (ali[1] != ali[0]).nonzero()[0]
            if len(muts):
                print adaID, fragment1+' - '+fragment2, ali.shape[1], muts









