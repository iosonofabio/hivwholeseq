# vim: fdm=marker
'''
author:     Fabio Zanini
date:       02/07/14
content:    Explore sites that are very conserved within our patient set, across
            the whole infection, compared to their behaviour within the subtype.
'''
# Modules
import os
import argparse
from itertools import izip, combinations, chain
from collections import defaultdict
from operator import itemgetter
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt
from Bio import SeqIO, AlignIO
from scipy.stats import pearsonr, spearmanr

from hivwholeseq.miseq import alpha, alphal
from hivwholeseq.patients.filenames import root_patient_folder
from hivwholeseq.sequencing.filenames import reference_folder
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.utils.sequence import translate_alignment
from hivwholeseq.utils.sequence import translate_with_gaps, get_subalignment
from hivwholeseq.one_site_statistics import (get_entropy,
    get_allele_frequencies_alignment)
from hivwholeseq.multipatient.get_shared_alleles_trajectories import (
    get_shared_allele_frequencies, get_patient_indices)
from hivwholeseq.analysis.explore_entropy_patsubtype import (
    get_subtype_reference_alignment, get_ali_entropy)



# Functions
def get_degenerate_dict():
    '''Get dictionary of degeneracies'''
    from collections import Counter
    from Bio.Data.CodonTable import standard_dna_table

    return dict(Counter(standard_dna_table.forward_table.itervalues()))


def get_codon_back_table():
    '''Get a complete back codon table'''
    from collections import defaultdict
    from Bio.Data.CodonTable import standard_dna_table
    
    table = defaultdict(list)
    for codon, aa in standard_dna_table.forward_table.iteritems():
        table[aa].append(codon)
    return dict(table)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Explore conservation levels across patients and subtype',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Regions to analyze (e.g. F1 p17)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot results')

    args = parser.parse_args()
    regions = args.regions
    VERBOSE = args.verbose
    plot = args.plot

    plt.ioff()
    data_all = {}

    for ir, region in enumerate(regions):
        if VERBOSE >= 1:
            print region

        if VERBOSE >= 2:
            print 'Get subtype B reference sequence alignment to HXB2 of', region
        ali_sub = get_subtype_reference_alignment(region, VERBOSE=VERBOSE)
        ali_subm = np.array(ali_sub, 'S1')

        if VERBOSE >= 2:
            print 'Get shared allele trajectories'
        # Low-coverage regions are bytecoded by -1
        data = get_shared_allele_frequencies(region, VERBOSE=VERBOSE, save=False)
        afts = np.ma.masked_less(data['afts'], -0.99)
        times = data['times']
        mapref = np.ma.masked_equal(data['mapref'], -1)
        refname = data['pnames'][-1]
        pnames = data['pnames'][:-1]
        samplenames = data['samplenames']
        n_patients = afts.shape[0]
        lseq = afts.shape[1]
        lalpha = afts.shape[2]
        patinds = get_patient_indices(samplenames, VERBOSE=VERBOSE)

        # Shrink coordinates (Pavel trimmed some protein e.g. RT)
        if mapref[-1] >= ali_sub.get_alignment_length():
            indpos = (mapref < ali_sub.get_alignment_length())
            mapref = mapref[indpos]
            afts = afts[:, :, indpos]
            lseq = afts.shape[1]


        if False:
            ##############################################
            # SITE ENTROPY
            ##############################################
            if VERBOSE >= 2:
                print 'Get entropy of patients'
            Spats = get_entropy(afts, VERBOSE=VERBOSE)


            if VERBOSE >= 2:
                print 'Get entropy in subtype alignment'
            Ssub = get_ali_entropy(ali_sub, positions=mapref, VERBOSE=VERBOSE)


            if VERBOSE >= 2:
                print 'Comparing entropy in each patient and in the subtype'
            data = defaultdict(list)
            Ssub_shuffled = Ssub.copy(); np.random.shuffle(Ssub_shuffled)


            # Time resolved
            for ipat, pname in enumerate(pnames):
                Spat = Spats[patinds[ipat]]
                data['spearmanr-time'].append([spearmanr(Spat_t, Ssub) for Spat_t in Spat])
                data['spearmanr-time-shuffled'].append([spearmanr(Spat_t, Ssub_shuffled)
                                                        for Spat_t in Spat])
            
            if plot:
                if VERBOSE >= 2:
                    print 'Plot entropy correlations with time'

                fig, ax = plt.subplots()
                for ipat, rps in enumerate(data['spearmanr-time']):
                    ts = times[patinds[ipat]]
                    ax.plot(ts, zip(*rps)[0], lw=2, label=pnames[ipat],
                            color=cm.jet(1.0 * ipat / len(pnames)))

                ax.set_xlabel('Time from infection [days]')
                ax.set_ylabel('Spearman r on entropy (pat VS subtype B)')

                ax.legend(loc=2, fontsize=10, ncol=2)
                ax.set_ylim(0, 0.5)
                ax.grid(True)
                ax.set_title(region)

                plt.tight_layout()
                plt.ion()
                plt.show()

        ###############################################
        ## DEGENERATE CODONS
        ###############################################
        #if VERBOSE >= 2:
        #    print 'Get translated alignment'
        #ali_prot = get_subtype_reference_alignment(region, VERBOSE=VERBOSE, type='aa')
        #if len(ali_prot) != len(ali_sub):
        #    # Sometimes the first sequence is HXB2 itself, for some obscure reason
        #    if ((len(ali_prot) == len(ali_sub) + 1) and 
        #        ('HXB2' in ali_prot[0].name) and \
        #        ([seq.name for seq in ali_sub] == [seq.name for seq in ali_prot[1:]])):
        #        ali_prot = ali_prot[1:]

        #    # If things go bad, just remake the translated alignment
        #    else:
        #        ali_prot = translate_alignment(ali_sub, VERBOSE=VERBOSE)

        ## Get 4-fold degenerate codon positions
        #pos_4fold = get_degenerate_pos()
        
        degenerate_dict = get_degenerate_dict()
        back_table = get_codon_back_table()

        data = defaultdict(list)
        data['entropy-4fold'] = [None for isample in xrange(afts.shape[0])]
        for isample, af in enumerate(afts):
            if VERBOSE >= 2:
                print isample + 1, 'of', afts.shape[0]

            # Skip incomplete for now (could be done better)
            if af[0].mask.any():
                continue

            consm = alpha[af.argmax(axis=0)]
            prot = translate_with_gaps(consm)
            prot_deg = np.array(map(degenerate_dict.get, prot), int)

            # NOTE: all 4-fold degenerate codons have all and only 3rd positions
            ind_4fold = (prot_deg == 4).nonzero()[0]
            ind_4fold_3rd = ind_4fold * 3 + 2
            Ssample = get_entropy(af[:, ind_4fold_3rd])

            # Get only subalignment that agrees at the first two positions
            Ssubsample = np.zeros(Ssample.shape[0])
            for ii, i in enumerate(ind_4fold):
                i *= 3
                #ind_sub = ((np.fromstring(ali_sub[:, i], 'S1') == consm[i]) &
                #           (np.fromstring(ali_sub[:, i + 1], 'S1') == consm[i + 1]))
                ind_sub = (ali_subm[:, i: i+2] == consm[i: i+2]).all(axis=1)
                ind_sub = ind_sub.nonzero()[0]
                ali_ssubm = ali_subm[ind_sub]
                nu3 = get_allele_frequencies_alignment(ali_ssubm[:, mapref[i+2]])
                Ssubsample[ii] = get_entropy(nu3, alphabet_axis=0)

            data['entropy-4fold'][isample] = (Ssample, Ssubsample)

        data['samplenames'] = samplenames
        data['times'] = times

        data_all[region] = dict(data)


    # Merge regions for calculating correlations
    data_merged = defaultdict(list)

    patients = load_patients()
    pnames = patients.index.tolist()

    samplenames_all = set.union(*[set(d['samplenames']) for d in data_all.itervalues()])
    fun = lambda x: pnames.index(x.split('_')[0]) * 100 + int(x.split('_')[1])
    samplenames_all = np.array(sorted(samplenames_all, key=fun))

    patinds = get_patient_indices(samplenames_all)

    data_merged['entropy-4fold'] = [None for i in xrange(len(samplenames_all))]
    data_merged['spearmanr'] = [None for i in xrange(len(samplenames_all))]
    for isample, samplename in enumerate(samplenames_all):
        Ssample = []
        for region, data in data_all.iteritems():
            if samplename not in data['samplenames']:
                break

            isn = data['samplenames'].tolist().index(samplename)
            Stmp = data['entropy-4fold'][isn]
            if Stmp is None:
                break

            Ssample.append(Stmp)

        else:
            Ssample = map(np.concatenate, izip(*Ssample))
            data_merged['entropy-4fold'][isample] = Ssample

            r = spearmanr(*Ssample)
            data_merged['spearmanr'][isample] = r

    if plot:
        if VERBOSE >= 2:
            print 'Plot 4fold degenerate entropy correlations with time'

        fig, ax = plt.subplots()
        for ipat, pname in enumerate(pnames):
            patient = Patient(patients.iloc[ipat])
            samplenames = samplenames_all[patinds[ipat]]
            indt = map(lambda x: int(x.split('_')[1]), samplenames)
            ts = patient.times[indt]

            rs = []
            for i in patinds[ipat]:
                tmp = data_merged['spearmanr'][i]
                if tmp is None:
                    rs.append(np.nan)
                else:
                    rs.append(tmp[0])

            ax.plot(ts, rs, lw=2, label=pnames[ipat],
                    color=cm.jet(1.0 * ipat / len(pnames)))

        ax.set_xlabel('Time from infection [days]')
        ax.set_ylabel('Spearman r on entropy (pat VS subtype B)')

        ax.legend(loc=2, fontsize=10, ncol=2)
        ax.set_ylim(0, 0.5)
        ax.grid(True)
        ax.set_title(' + '.join(regions))

        plt.tight_layout()
        plt.ion()
        plt.show()

