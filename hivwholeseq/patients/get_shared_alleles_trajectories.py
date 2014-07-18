# vim: fdm=marker
'''
author:     Fabio Zanini
date:       02/07/14
content:    Get trajectories of alleles across patients.
'''
# Modules
import os
import argparse
from operator import itemgetter
import numpy as np
from Bio import SeqIO, AlignIO

from hivwholeseq.miseq import alpha
from hivwholeseq.mapping_utils import align_muscle
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.filenames import get_initial_consensus_filename, \
        get_mapped_to_initial_filename, get_allele_frequency_trajectories_filename, \
        get_allele_count_trajectories_filename
from hivwholeseq.patients.one_site_statistics import \
        plot_allele_frequency_trajectories_from_counts as plot_nus_from_act
from hivwholeseq.patients.one_site_statistics import \
        plot_allele_frequency_trajectories_from_counts_3d as plot_nus_from_act_3d
from hivwholeseq.patients.one_site_statistics import get_allele_count_trajectories
from hivwholeseq.fork_cluster import fork_get_allele_frequency_trajectory as fork_self



# Function
def build_coordinate_maps(ali, VERBOSE=0, stripgaps=False):
    '''Build coordinate maps from alignment to single references'''

    rows = len(ali)
    cols = ali.get_alignment_length()

    maps = np.ma.masked_all((rows, cols), int)

    for i, seq in enumerate(ali):
        smat = np.array(seq, 'S1')
        is_nongap = smat != '-'
        maps[i, is_nongap] = np.arange((is_nongap).sum())

    if stripgaps:
        maps = maps[:, -(maps.mask.any(axis=0))]

    return maps



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get shared allele trajectories',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the allele frequency trajectories')
    parser.add_argument('--logit', action='store_true',
                        help='use logit scale (log(x/(1-x)) in the plots')
    parser.add_argument('--PCR1', type=int, default=1,
                        help='Take only PCR1 samples [0=off, 1=where both available, 2=always]')

    args = parser.parse_args()
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    plot = args.plot
    use_PCR1 = args.PCR1
    use_logit = args.logit


    patients = load_patients()
    # FIXME: extend to all patients
    patients = patients.loc[-patients.index.isin(['15363', '15376', '15313', '20529', '15034', '20097'])]
    if pnames is not None:
        patients = patients.loc[patients.index.isin(pnames)]

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    for fragment in fragments:
        if VERBOSE >= 1:
            print fragment

        if VERBOSE >= 1:
            print 'Collect initial references'
        refs = []
        for pname, patient in patients.iterrows():
            patient = Patient(patient)
            ref = patient.get_reference(fragment)
            refs.append(ref)

        if VERBOSE >= 1:
            print 'Align references'
        ali = align_muscle(*refs, sort=True)

        if VERBOSE >= 1:
            print 'Getting coordinate maps'
        # Exclude sites that are not present in all patients
        maps = build_coordinate_maps(ali, VERBOSE=VERBOSE, stripgaps=True)

        afts = np.zeros((maps.shape[0], maps.shape[1], len(alpha)), object)

        if VERBOSE >= 1:
            print 'Collecting alleles'
        for ip, (pname, patient) in enumerate(patients.iterrows()):
            patient = Patient(patient)
            patient.discard_nonsequenced_samples()
            samplenames = patient.samples.index

            # Collect allele counts from patient samples, and return only positive hits
            # sns contains sample names and PCR types
            (sns, act) = get_allele_count_trajectories(pname, samplenames, fragment,
                                                       use_PCR1=use_PCR1, VERBOSE=VERBOSE)
            ind = [i for i, (_, sample) in enumerate(patient.samples.iterrows())
                   if sample.name in map(itemgetter(0), sns)]
            samples = patient.samples.iloc[ind]
            times = (samples.date - patient.transmission_date) / np.timedelta64(1, 'D')
            # FIXME: determine better PCR errors?
            ntemplate = sample['n templates']
            if ntemplate > 0.1:
                depthmax = np.maximum(1e-3, 1.0 / ntemplate)
            else:
                depthmax = 1e-3

            # Low-coverage sampled are bitcoded as -1
            aft = np.zeros_like(act, dtype=float)
            for i in xrange(aft.shape[0]):
                for k in xrange(aft.shape[2]):
                    ac = act[i, :, k]
                    co = ac.sum()
                    if co < 1000:
                        aft[i, :, k] = -1
                        continue

                    af = 1.0 * ac / co
                    af[af < depthmax] = 0
                    af[af > 1 - depthmax] = 1
                    af /= af.sum()
                    aft[i, :, k] = af

            mapi = maps[ip]
            for i, ind in enumerate(mapi):
                for j in xrange(len(alpha)):
                    afts[ip, i, j] = aft[:, j, ind]

        from hivwholeseq.patients.filenames import root_patient_folder
        afts.dump(root_patient_folder+'all/aft_shared_'+fragment+'.npy')
        AlignIO.write(ali, root_patient_folder+'all/aft_shared_ali_'+fragment+'.fasta', 'fasta')
        maps.dump(root_patient_folder+'all/aft_shared_maps_'+fragment+'.npy')
