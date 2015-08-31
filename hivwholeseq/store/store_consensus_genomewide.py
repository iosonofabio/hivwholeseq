# vim: fdm=marker
'''
author:     Fabio Zanini
date:       24/10/13
content:    Build genome wide reference from the 6 fragments.
'''
# Modules
import os
import sys
import argparse
from itertools import izip
import numpy as np
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO as SeqIO

from hivwholeseq.patients.patients import SamplePat
from hivwholeseq.patients.patients import load_samples_sequenced as lssp
from hivwholeseq.sequencing.check_overlaps import get_overlapping_fragments



# Functions
def merge_fragments(sequences, name='', VERBOSE=0):
    '''Merge references at overlapping pairs'''
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.Alphabet.IUPAC import ambiguous_dna
    from seqanpy import align_ladder

    from hivwholeseq.utils.sequence import pretty_print_pairwise_ali

    consensus = []
    seq_old = ''.join(sequences['F1'])
    for i in xrange(5):
        seq_new = ''.join(sequences['F'+str(i+2)])
        (score, ali1, ali2) = align_ladder(seq_old, seq_new, score_gapopen=-10)

        if VERBOSE >= 3:
            pretty_print_pairwise_ali([ali1, ali2], name1='F'+str(i+1), name2='F'+str(i+2))

        # Overlap: the first sequence is better at the start, the second at the end
        end1 = len(ali1.rstrip('-'))
        start2 = len(ali2) - len(ali2.lstrip('-'))
        len_overlap = end1 - start2

        # There might a too short consensus, just join them with N
        if len_overlap < 50:
            consensus.append(seq_old)
            consensus.append('N' * 10)
            if i == 4:
                consensus.append(seq_new)
            else:
                seq_old = seq_new

            continue

        overlap1 = np.fromstring(ali1[start2: end1], 'S1')
        overlap2 = np.fromstring(ali2[start2: end1], 'S1')
        overlap = overlap1.copy()
        ind_overlap_mismatch = (overlap1 != overlap2).nonzero()[0]
        for j in ind_overlap_mismatch:
            if j < len(overlap) // 3:
                continue
            elif j < 2 * len(overlap) // 3:
                overlap[j] = 'N'
            else:
                overlap[j] = overlap2[j]
        overlap = overlap.tostring()

        consensus.append(ali1[:start2])
        consensus.append(overlap)
        if i == 4:
            consensus.append(ali2[end1:])
        else:
            seq_old = ali2[end1:].replace('-', '')


    consensus = ''.join(consensus)
    cons_rec = SeqRecord(Seq(consensus, IUPAC.ambiguous_dna),
                         id=name, name=name,
                         description=name+', genomewide')

    return cons_rec



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Merge consensi of the 6 fragments')
    pats_or_samples = parser.add_mutually_exclusive_group(required=True)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to map')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--save', action='store_true',
                        help='Save genomewide consensi to file')
    
    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    VERBOSE = args.verbose
    use_save = args.save

    samples_pat = lssp()
    if pnames is not None:
        samples_pat = samples_pat.loc[samples_pat.patient.isin(pnames)]
    elif samplenames is not None:
        samples_pat = samples_pat.loc[samples_pat.index.isin(samplenames)]

    if VERBOSE >= 2:
        print 'samples', samples_pat.index.tolist()

    consensi = []
    for samplename_pat, sample_pat in samples_pat.iterrows():
        sample_pat = SamplePat(sample_pat)
        if VERBOSE >= 1:
            print sample_pat.name,

        consensi_pat = {}
        try:
            for i in xrange(6):
                fragment = 'F'+str(i+1)
                consensi_pat[fragment] = sample_pat.get_consensus(fragment)

        except IOError:
            if VERBOSE >= 1:
                print 'warning: some consensus not found: skipping'
            continue

        except ValueError:
            if VERBOSE >= 1:
                print 'warning: some consensus is not clean (e.g. contamination): skipping'
            continue

        if VERBOSE >= 1:
            print 'ok: all 6 consensi found'

        consensus =  merge_fragments(consensi_pat, name=sample_pat.patient,
                                     VERBOSE=VERBOSE)

        if len(consensus) < 8800:
            print 'WARNING: the consensus looks too short!'

        if use_save:
            output_filename = sample_pat.get_consensus_filename('genomewide')
            SeqIO.write(consensus, output_filename, 'fasta')
            if VERBOSE >= 1:
                print 'Genomewide consensus written'

