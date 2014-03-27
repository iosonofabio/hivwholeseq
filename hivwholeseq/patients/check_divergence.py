# vim: fdm=marker
'''
author:     Fabio Zanini
date:       26/03/14
content:    Check divergence of sequences during time in various ways.
'''
# Modules
import os
import argparse
import numpy as np

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.patients import get_patient
from hivwholeseq.patients.filenames import get_initial_consensus_filename, \
        get_mapped_to_initial_filename, get_allele_frequency_trajectories_filename, \
        get_allele_count_trajectories_filename



# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Check patient samples')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose

    patient = get_patient(pname)
    times = patient.times()

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    if VERBOSE:
        print pname
    divs = []
    for fragment in fragments:

        # 1: use allele frequency trajectories
        aft_filename = get_allele_frequency_trajectories_filename(pname, fragment)
        aft = np.ma.array(np.load(aft_filename))
        aft[aft < 1e-2] = np.ma.masked

        cons = alpha[aft[0].argmax(axis=0)]
        anc = frozenset(zip(cons, np.arange(len(cons))))

        from itertools import izip
        n_muts = []
        for af in aft:
            muts = frozenset([(alpha[ai], pos)
                              for (ai, pos) in izip(*(af > 0.5).nonzero())
                              if alpha[ai] != '-'])
            muts -= anc
            n_muts.append(len(muts))

        div = (0.5 * n_muts[-1] / times[-1] / len(cons))
        div += (0.5 * n_muts[-2] / times[-2] / len(cons))
        divs.append(div)
        if VERBOSE:
            print fragment, '{:1.2e}'.format(div * 365.25), 'substitutions per site per year'

        # Divergence at low frequency
        nutots = []
        for af in aft:
            nutot = sum(1 for (ai, pos) in izip(*((af > 0.01) & (af < 0.1)).nonzero())
                        if (alpha[ai] != '-') and (cons[pos] != alpha[ai]))
            nutots.append(nutot)
        print fragment, 'nutots', nutots 


