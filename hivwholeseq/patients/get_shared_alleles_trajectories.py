# vim: fdm=marker
'''
author:     Fabio Zanini
date:       02/07/14
content:    Get trajectories of alleles across patients.
'''
# Modules
import os
import argparse
from collections import defaultdict
from operator import itemgetter
import numpy as np
from Bio import SeqIO, AlignIO

from hivwholeseq.miseq import alpha
from hivwholeseq.mapping_utils import align_muscle
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_mapped_to_initial_filename, get_allele_frequency_trajectories_filename, \
        get_allele_count_trajectories_filename
from hivwholeseq.patients.one_site_statistics import \
        plot_allele_frequency_trajectories_from_counts as plot_nus_from_act
from hivwholeseq.patients.one_site_statistics import \
        plot_allele_frequency_trajectories_from_counts_3d as plot_nus_from_act_3d
from hivwholeseq.patients.one_site_statistics import get_allele_count_trajectories
from hivwholeseq.cluster.fork_cluster import fork_get_allele_frequency_trajectory as fork_self



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


def build_coordinate_map_reference(alimap, refg, VERBOSE=0):
    '''Figure out where the aligned alleles are in the reference, for convenience'''
    refg = np.array(refg)
    refcs = (refg != '-').cumsum() - 1
    refcoo = np.array([refcs[pos] for pos in alimap], int)
    return refcoo


def get_shared_allele_frequencies(region, pnames=None, VERBOSE=0, save=False,
                                  reference='HXB2'):
    '''Align allele frequencies from several patients in a region'''
    import os
    from hivwholeseq.patients.filenames import root_patient_folder

    fn_aft = root_patient_folder+'all/aft_shared_'+region+'.npz'
    fn_ali = root_patient_folder+'all/aft_shared_ali_'+region+'.fasta'
    fn_map = root_patient_folder+'all/aft_shared_maps_'+region+'.npz'

    # Recycle existing files if possible
    if (not save) and all(map(os.path.isfile, (fn_aft, fn_ali, fn_map))):
        npdata = np.load(fn_aft)
        ali = AlignIO.read(fn_ali, 'fasta')
        mapdict = np.load(fn_map)

        data = dict(npdata)
        data['ali'] = ali
        data.update(mapdict)

    else:
        if VERBOSE >= 1:
            print 'Load patients'
        from hivwholeseq.patients.patients import load_patients, Patient
        patients = load_patients()
        if pnames is not None:
            patients = patients.loc[patients.index.isin(pnames)]

        if VERBOSE >= 1:
            print 'Collect initial references'
        refs = []
        for pname, patient in patients.iterrows():
            patient = Patient(patient)
            ref = patient.get_reference(region)
            refs.append(ref)
        if reference is not None:
            from hivwholeseq.reference import load_custom_reference
            refs.append(load_custom_reference(reference, region=region))

        if VERBOSE >= 1:
            print 'Align references'
        ali = align_muscle(*refs, sort=True)

        if VERBOSE >= 1:
            print 'Getting coordinate maps'
        # Exclude sites that are not present in all patients + reference
        # NOTE: the reference does not really restrict and is very convenient
        maps = build_coordinate_maps(ali, VERBOSE=VERBOSE, stripgaps=True)

        # Get map to reference instead of alignment, for convenience
        if reference is not None:
            mapref = build_coordinate_map_reference(maps[-1], ali[-1], VERBOSE=VERBOSE)

        if VERBOSE >= 1:
            print 'Collecting alleles'
        data = defaultdict(list)
        for ip, (pname, patient) in enumerate(patients.iterrows()):
            patient = Patient(patient)
            patient.discard_nonsequenced_samples()

            # Collect allele counts from patient samples
            act, indt = patient.get_allele_count_trajectories(region,
                                                             VERBOSE=VERBOSE,
                                                             safe=True)
            timespat = patient.times[indt]
            depthmax = np.maximum(1.0 / patient.n_templates[indt], 2e-3)
            #FIXME: for now, put 2e-3 to the masked positions, but this is no good
            depthmax[depthmax.mask] = 2e-3

            # Low-coverage sampled are bytecoded as -1
            aft = np.zeros_like(act, dtype=float)
            for i in xrange(aft.shape[0]):
                for k in xrange(aft.shape[2]):
                    ac = act[i, :, k]
                    co = ac.sum()
                    if co < 1000:
                        aft[i, :, k] = -1
                        continue

                    af = 1.0 * ac / co
                    af[af < depthmax[i]] = 0
                    af[af > 1 - depthmax[i]] = 1
                    af /= af.sum()
                    aft[i, :, k] = af

            # Set the data
            data['depthmaxs'].append(depthmax.data)
            data['times'].append(timespat)
            mapi = maps[ip]
            afts = np.zeros((aft.shape[0], len(alpha), len(mapi)), float)
            for i, ind in enumerate(mapi):
                afts[:, :, i] = aft[:, :, ind]
            data['afts'].append(afts)
            data['samplenames'].append([pname+'_'+str(it) for it in indt])

        # Stack
        data = dict(data)
        for key, val in data.iteritems():
            data[key] = np.concatenate(val, axis=0)

        # Idem bytecode for maps
        maps_bytecode = np.ma.filled(maps, -1)

        mapdict = {'maps': maps_bytecode,
                   'pnames': np.append(np.array(patients.index), reference)}
        if reference is not None:
            mapdict['mapref'] = mapref

        if save:
            if VERBOSE >= 1:
                print 'Saving to file'
            np.savez_compressed(fn_aft, **data)
            AlignIO.write(ali, fn_ali, 'fasta')
            np.savez_compressed(fn_map, **mapdict)

        data.update(mapdict)
        data['ali'] = ali

    return data


def get_patient_indices(samplenames, VERBOSE=0):
    '''Get which sample is in which patient'''
    inds = []
    ind = None
    pname_old = None
    for i, samplename in enumerate(samplenames):
        pname = samplename.split('_')[0]
        if pname == pname_old:
            ind.append(i)
        else:
            inds.append(ind)
            ind = [i]
            pname_old = pname
    inds.append(ind)
    inds = inds[1:]

    return inds



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get shared allele trajectories',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--regions', nargs='*',
                        help='Regions to analyze (e.g. F1 V3)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--reference', default='HXB2',
                        help='External reference to align')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    reference = args.reference

    if not regions:
        regions = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'regions', regions

    for region in regions:
        if VERBOSE >= 1:
            print region

        data = get_shared_allele_frequencies(region, pnames, VERBOSE=VERBOSE,
                                             reference=reference,
                                             save=True)


