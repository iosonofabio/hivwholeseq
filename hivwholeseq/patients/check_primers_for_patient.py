# vim: fdm=marker
'''
author:     Fabio Zanini
date:       24/05/14
content:    Use the overlap to check on conservation of primers, to decide for
            future samples.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import reverse_complement as rc

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.miseq import alpha
from hivwholeseq.patients.filenames import get_initial_consensus_filename, \
        get_allele_count_trajectories_filename
from hivwholeseq.patients.patients import get_patient
from hivwholeseq.primer_info import primers_outer, find_primer_seq
from hivwholeseq.sequence_utils import expand_ambiguous_seq



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Update initial consensus')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='+',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose

    patient = get_patient(pname)

     # If the script is called with no fragment, iterate over all
    # NOTE: This assumes that the fragment has been sequenced. If not, we should
    # take a later time point or similia, but let's not worry about it yet.
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments
   
    
    for fragment in fragments:
        print fragment
        # Check the overlap from the left for fwd primer, unless F1
        if fragment != 'F1':
            frag_left = 'F'+str(int(fragment[1]) - 1)
            cons_fn_left = get_initial_consensus_filename(pname, frag_left)
            cons_left = SeqIO.read(cons_fn_left, 'fasta')

            if fragment == 'F5':
                frag_spec = 'F5a'
            else:
                frag_spec = fragment

            primer_fwd = primers_outer[frag_spec][0]
            primer_fwd_exps = expand_ambiguous_seq(primer_fwd)

            pos = find_primer_seq(cons_left, primer_fwd)
            primer_fwd_pat = str(cons_left[pos: pos + len(primer_fwd)].seq)

            if primer_fwd_pat in primer_fwd_exps:
                status = 'OK'
            else:
                status = 'FAIL'
                primer_fwd_expsm = np.array([np.fromstring(pre, 'S1') for pre in primer_fwd_exps], 'S1')
                primer_fwd_closest = primer_fwd_exps[(primer_fwd_expsm == primer_fwd_pat).argmin()]

            print 'FWD:', status
            print 'Designed:\t\t', primer_fwd
            if status == 'FAIL':
                print 'Design (closest)\t', primer_fwd_closest
            print 'Patient:\t\t', primer_fwd_pat

            if status == 'FAIL':
                act_filename = get_allele_count_trajectories_filename(pname, frag_left)
                act = np.load(act_filename)
                ac_pr = act[0, :, pos: pos + len(primer_fwd)]
                ind = (np.fromstring(primer_fwd_closest, 'S1') != np.fromstring(primer_fwd_pat, 'S1')).nonzero()[0]
                for i in ind:
                    print i, primer_fwd_closest[i], primer_fwd_pat[i], ac_pr[:, i]

        # Check overlap from the right for rev primer
        if fragment != 'F6':
            frag_right = 'F'+str(int(fragment[1]) + 1)
            cons_fn_right = get_initial_consensus_filename(pname, frag_right)
            cons_right = SeqIO.read(cons_fn_right, 'fasta')

            if fragment == 'F5':
                frag_spec = 'F5a'
            else:
                frag_spec = fragment

            primer_rev_rc= primers_outer[frag_spec][1]
            primer_rev = rc(primer_rev_rc)
            primer_rev_exps = expand_ambiguous_seq(primer_rev)

            pos = find_primer_seq(cons_right, primer_rev_rc)
            primer_rev_pat = rc(str(cons_right[pos: pos + len(primer_rev_rc)].seq))

            if primer_rev_pat in primer_rev_exps:
                status = 'OK'
            else:
                status = 'FAIL'
                primer_rev_expsm = np.array([np.fromstring(pre, 'S1') for pre in primer_rev_exps], 'S1')
                primer_rev_closest = primer_rev_exps[(primer_rev_expsm == primer_rev_pat).argmin()]

            print 'REV:', status
            print 'Designed:\t\t', primer_rev
            if status == 'FAIL':
                print 'Design (closest):\t', primer_rev_closest
            print 'Patient:\t\t', primer_rev_pat

            if status == 'FAIL':
                act_filename = get_allele_count_trajectories_filename(pname, frag_right)
                act = np.load(act_filename)
                ac_pr = act[0, :, pos: pos + len(primer_rev)]
                tmp = ac_pr[0].copy()
                ac_pr[0] = ac_pr[3]
                ac_pr[3] = tmp
                tmp = ac_pr[1].copy()
                ac_pr[1] = ac_pr[2]
                ac_pr[2] = tmp
                ac_pr = ac_pr[:, ::-1]

                ind = (np.fromstring(primer_rev_closest, 'S1') != np.fromstring(primer_rev_pat, 'S1')).nonzero()[0]
                for i in ind:
                    print i, primer_rev_closest[i], primer_rev_pat[i], ac_pr[:, i]

        print ''
