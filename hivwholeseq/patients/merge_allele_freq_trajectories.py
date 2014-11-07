# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/03/14
content:    Merge allele frequency trajectories of all fragments.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO

from seqanpy import align_overlap

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_mapped_to_initial_filename, get_allele_frequency_trajectories_filename, \
        get_allele_count_trajectories_filename
from hivwholeseq.patients.one_site_statistics import \
        plot_allele_frequency_trajectories_from_counts as plot_nus_from_act
from hivwholeseq.patients.one_site_statistics import \
        plot_allele_frequency_trajectories_from_counts_3d as plot_nus_from_act_3d



# Functions
def merge_allele_count_trajectories(ref_genomewide, acss, VERBOSE=0):
    '''Merge the allele frequency trajectories of all fragments by summing counts'''

    _, _, inds = zip(*acss)
    ind_min = min(map(min, inds))
    ind_max = max(map(max, inds))
    acs = np.zeros((ind_max + 1 - ind_min, len(alpha), len(ref_genomewide)), int)

    pos_ref = 1000
    for ifr, (ref, acsi, ind) in enumerate(acss):
        ind -= ind_min

        # Find the coordinates
        (score, ali1, ali2) = align_overlap(ref_genomewide[pos_ref - 1000: pos_ref + 3000],
                                            ref, score_gapopen=-20)
        fr_start = len(ali2) - len(ali2.lstrip('-'))
        fr_end = len(ali2.rstrip('-'))

        if VERBOSE:
            print 'F'+str(ifr+1), pos_ref - 1000 + fr_start, pos_ref - 1000 + fr_end

        # Scan the alignment
        pos_ref = pos_ref - 1000 + fr_start
        pos_fr = 0
        for pos_ali in xrange(fr_start, fr_end):
            # Gap in genomewise, ignore position
            if ali1[pos_ali] == '-':
                pos_fr += 1
                continue

            # Gap in fragment, ignore FIXME: probably we should put deletions! 
            elif ali2[pos_ali] == '-':
                pos_ref += 1
                continue

            acs[ind, :, pos_ref] = acsi[:, :, pos_fr]
            pos_fr += 1
            pos_ref += 1

    ind = np.arange(ind_min, ind_max + 1)

    return (acs, ind)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Update initial consensus',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--save', action='store_true',
                        help='Save the allele frequency trajectories to file')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the allele frequency trajectories')
    parser.add_argument('--logit', action='store_true',
                        help='use logit scale (log(x/(1-x)) in the plots')
    parser.add_argument('--PCR1', type=int, default=1,
                        help='Take only PCR1 samples [0=off, 1=where both available, 2=always]')

    args = parser.parse_args()
    pnames = args.patients
    VERBOSE = args.verbose
    save_to_file = args.save
    plot = args.plot
    use_PCR1 = args.PCR1
    use_logit = args.logit

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    fragments = ['F'+str(i) for i in xrange(1, 7)]
    
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        if VERBOSE:
            print pname

        # Get the initial genomewide consensus
        conss_genomewide = ''.join(patient.get_reference('genomewide'))

        # Collect the allele count trajectories
        acts = []
        for fragment in fragments:
            ref = ''.join(patient.get_reference(fragment))
            act, ind = patient.get_allele_count_trajectories(fragment)
            acts.append((ref, act, ind))

        # Merge allele counts
        (act, ind) = merge_allele_count_trajectories(conss_genomewide, acts,
                                                     VERBOSE=VERBOSE)
        # Normalize to frequencies
        afs = (1.0 * act.swapaxes(0, 1) / (0.1 + act.sum(axis=1))).swapaxes(0, 1)

        if save_to_file:
            fn_out = get_allele_count_trajectories_filename(pname, 'genomewide')
            np.savez(fn_out, ind=ind, act=act)

        if plot is not None:
            import matplotlib.pyplot as plt

            times = patient.times[ind]
            ntemplates = patient.n_templates[ind]

            if plot in ('2D', '2d', ''):
                plot_nus_from_act(times, act, title='Patient '+pname, VERBOSE=VERBOSE,
                                  ntemplates=ntemplates,
                                  logit=use_logit,
                                  threshold=0.9)

            if plot in ('3D', '3d'):
                plot_nus_from_act_3d(times, act, title='Patient '+pname, VERBOSE=VERBOSE,
                                     logit=use_logit,
                                     threshold=0.9)

    if plot is not None:
        plt.ion()
        plt.show()
