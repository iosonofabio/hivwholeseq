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
from itertools import izip
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio import SeqIO
from Bio.Seq import reverse_complement as rc

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_allele_count_trajectories_filename, get_primers_filename
from hivwholeseq.patients.patients import load_patient, SamplePat
from hivwholeseq.patients.samples import load_sample_sequenced
from hivwholeseq.sequencing.primer_info import primers_outer, find_primer_seq
from hivwholeseq.sequence_utils import expand_ambiguous_seq



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Update initial consensus')
    pat_or_sample = parser.add_mutually_exclusive_group(required=True)
    pat_or_sample.add_argument('--patient',
                               help='Patient to analyze')
    pat_or_sample.add_argument('--sample',
                               help='Sample to analyze')
    parser.add_argument('--fragments', nargs='+',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--save', action='store_true',
                        help='Save patient primers to file')

    args = parser.parse_args()
    pname = args.patient
    samplename = args.sample
    fragments = args.fragments
    VERBOSE = args.verbose
    use_save = args.save

    if pname is not None:
        patient = load_patient(pname)
        seq_fun = patient.get_reference
    elif samplename is not None:
        sample = load_sample_sequenced(samplename)
        seq_fun = sample.get_consensus

     # If the script is called with no fragment, iterate over all
    # NOTE: This assumes that the fragment has been sequenced. If not, we should
    # take a later time point or similia, but let's not worry about it yet.
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments
    
    primers_std = defaultdict(dict)
    primers_clo = defaultdict(dict)
    primers_pat = defaultdict(dict)
    for fragment in fragments:
        if pname is not None:
            print pname, fragment
        elif samplename is not None:
            print samplename, fragment

        # Check the overlap from the left for fwd primer, unless F1
        if fragment != 'F1':
            frag_left = 'F'+str(int(fragment[1]) - 1)
            cons_left = seq_fun(frag_left)

            # Only use the best-performing primers...
            if fragment == 'F5':
                frag_spec = 'F5a'
            elif fragment == 'F3':
                frag_spec = 'F3b'
            else:
                frag_spec = fragment

            primer_fwd = primers_outer[frag_spec][0]
            primer_fwd_exps = expand_ambiguous_seq(primer_fwd)

            pos = find_primer_seq(cons_left, primer_fwd)
            primer_fwd_pat = str(cons_left[pos: pos + len(primer_fwd)].seq)
            primer_fwd_patm = np.fromstring(primer_fwd_pat, 'S1')

            primer_fwd_expsm = np.array([np.fromstring(pre, 'S1') for pre in primer_fwd_exps], 'S1')
            primer_fwd_closest = primer_fwd_exps[(primer_fwd_expsm == primer_fwd_patm).sum(axis=1).argmax()]

            if primer_fwd_pat in primer_fwd_exps:
                status = 'OK'
            else:
                status = 'FAIL'

            print 'FWD:', status
            print 'Designed:\t\t', '5\'', primer_fwd, '3\''
            if status == 'FAIL':
                print 'Design (closest):\t', '5\'', primer_fwd_closest, '3\''
                print '                 \t  ', 
                match_string = []
                for nuc_des, nuc_pat in izip(primer_fwd_closest, primer_fwd_pat):
                    if nuc_des != nuc_pat:
                        match_string.append('x')
                    else:
                        match_string.append(' ')
                print ''.join(match_string)

            print 'Patient:\t\t', '5\'', primer_fwd_pat, '3\''

            primers_std[fragment]['fwd'] = primer_fwd
            primers_clo[fragment]['fwd'] = primer_fwd_closest
            primers_pat[fragment]['fwd'] = primer_fwd_pat

            if (pname is not None) and (status == 'FAIL'):
                sample = SamplePat(patient.samples.iloc[0])
                ac_filename = sample.get_allele_counts_filename(frag_left)
                if not os.path.isfile(ac_filename):
                    print 'Allele count filename not found.'
                else:
                    ac = np.load(ac_filename).sum(axis=0)
                    ac_pr = ac[:, pos: pos + len(primer_fwd)]
                    ind = (np.fromstring(primer_fwd_closest, 'S1') != np.fromstring(primer_fwd_pat, 'S1')).nonzero()[0]
                    for i in ind:
                        print i, primer_fwd_closest[i], primer_fwd_pat[i], ac_pr[:, i]

        # Check overlap from the right for rev primer
        if fragment != 'F6':
            frag_right = 'F'+str(int(fragment[1]) + 1)
            cons_right = seq_fun(frag_right)

            if fragment == 'F5':
                frag_spec = 'F5a'
            elif fragment == 'F3':
                frag_spec = 'F3b'
            else:
                frag_spec = fragment

            primer_rev_rc= primers_outer[frag_spec][1]
            primer_rev = rc(primer_rev_rc)
            primer_rev_exps = expand_ambiguous_seq(primer_rev)

            pos = find_primer_seq(cons_right, primer_rev_rc)
            primer_rev_pat = rc(str(cons_right[pos: pos + len(primer_rev_rc)].seq))
            primer_rev_patm = np.fromstring(primer_rev_pat, 'S1')

            primer_rev_expsm = np.array([np.fromstring(pre, 'S1') for pre in primer_rev_exps], 'S1')
            primer_rev_closest = primer_rev_exps[(primer_rev_expsm == primer_rev_patm).sum(axis=1).argmax()]

            if primer_rev_pat in primer_rev_exps:
                status = 'OK'
            else:
                status = 'FAIL'

            print 'REV:', status
            print 'Designed:\t\t', '5\'', primer_rev, '3\''
            if status == 'FAIL':
                print 'Design (closest):\t', '5\'', primer_rev_closest, '3\''
                print '                 \t  ', 
                match_string = []
                for nuc_des, nuc_pat in izip(primer_rev_closest, primer_rev_pat):
                    if nuc_des != nuc_pat:
                        match_string.append('x')
                    else:
                        match_string.append(' ')
                print ''.join(match_string)
            print 'Patient:\t\t', '5\'', primer_rev_pat, '3\''

            primers_std[fragment]['rev'] = primer_rev
            primers_clo[fragment]['rev'] = primer_rev_closest
            primers_pat[fragment]['rev'] = primer_rev_pat

            if (pname is not None) and (status == 'FAIL'):
                sample = SamplePat(patient.samples.iloc[0])
                ac_filename = sample.get_allele_counts_filename(frag_right)
                if not os.path.isfile(ac_filename):
                    print 'Allele count filename not found.'
                else:
                    ac = np.load(ac_filename).sum(axis=0)
                    ac_pr = ac[:, pos: pos + len(primer_rev)]
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

    if use_save and (pname is not None):
        fn = get_primers_filename(pname, format='fasta')
        
        seqs = []
        for fragment in primers_pat:
            for sense in primers_pat[fragment]:
                pr_std = primers_std[fragment][sense]
                pr_clo = primers_clo[fragment][sense]
                pr_pat = primers_pat[fragment][sense]
                
                pr_std = SeqRecord(Seq(pr_std, ambiguous_dna),
                                   id=fragment+sense+'StandardPrimer',
                                   name=fragment+sense+'StandardPrimer',
                                   description=fragment+', '+sense+', standard primer',
                                  )
                seqs.append(pr_std)

                pr_clo = SeqRecord(Seq(pr_clo, ambiguous_dna),
                                   id=fragment+sense+'StandardPrimerClosest',
                                   name=fragment+sense+'StandardPrimerClosest',
                                   description=fragment+', '+sense+', standard primer, closest non-ambiguous',
                                  )
                seqs.append(pr_clo)

                pr_pat = SeqRecord(Seq(pr_pat, ambiguous_dna),
                                   id=fragment+sense+'Pat'+pname,
                                   name=fragment+sense+'Pat'+pname,
                                   description=fragment+', '+sense+', patient '+pname,
                                  )
                seqs.append(pr_pat)

        SeqIO.write(seqs, fn, 'fasta')


