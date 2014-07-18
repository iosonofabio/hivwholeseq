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
from hivwholeseq.filenames import reference_folder
from hivwholeseq.patients.patients import load_patients, Patient

# Globals
tree_ali_foldername = reference_folder+'alignments_trees/typeM/'
colors = {'p17': 'r',
          'p24': 'g',
          'PR': 'r',
          'RT': 'g',
          'p15': 'purple',
          'IN': 'orange',
          'gp41': 'r'}
protnames = {'F1': ['p17', 'p24'],
             'F2': ['PR', 'RT'],
             'F3': ['p15', 'IN'],
             'F4': [],
             'F5': [],
             'F6': ['gp41']}
protnames_sort = list(chain(*[protnames['F'+str(i+1)] for i in xrange(6)]))



# Functions
def get_subtype_B_obs(protname, seqm, VERBOSE=0):
    '''Get the observables from subtype B alignments'''
    aliB_fn = tree_ali_foldername+protname+'.B.nuc.aligned.fasta'
    aliB = AlignIO.read(aliB_fn, 'fasta')
    aliBm = np.array(aliB, 'S1')
    # Get allele freqs ignoring ambiguous: N, Y, R, etc. and renormalizing
    nuBm = np.array([(aliBm == a).mean(axis=0) for a in alpha[:5]])
    nuBm /= nuBm.sum(axis=0)
    seqBm = alpha[nuBm.argmax(axis=0)]
    seqBm[seqBm == '-'] = 'X'
    S = np.zeros(len(seqBm))
    for nuB in nuBm:
        S -= (nuB + 1e-8) * np.log2(nuB + 1e-8)

    # Align alignments -- pairwise ;-)
    band = 50

    # Start
    seed = seqBm[:10]
    n_match = [(seqm[i: i + len(seed)] == seed).sum() for i in xrange(len(seqm) - len(seed))]
    pos = np.argmax(n_match)
    if n_match[pos] < 0.65 * len(seed):
        raise ValueError('ERROR: start not found')
    start = pos

    # End
    if protname == 'gp41':
        threshold = 0.6
    else:
        threshold = 0.7
    seed = seqBm[-10:]
    lp = len(seqBm)
    search_start = start + lp - band
    search_end = min(search_start + 2 * band, len(seqm) - len(seed))
    n_match = [(seqm[i: i + len(seed)] == seed).sum() for i in xrange(search_start, search_end)]
    pos = np.argmax(n_match)
    if n_match[pos] < threshold * len(seed):
        # Integrase ends slightly after the end of F3a, so we try and do the opposite, we trim IN
        if protname == 'IN':
            end = len(seqm)
            seed = seqm[-15:]
            search_start = len(seqBm) - 100
            search_end = len(seqBm) - len(seed)
            n_match = [(seqBm[i: i + len(seed)] == seed).sum() for i in xrange(search_start, search_end)]
            pos = np.argmax(n_match)
            if n_match[pos] < 0.7 * len(seed):
                raise ValueError('ERROR: neither end of IN in F3, nor end of F3 in IN found')
            endB = search_start + pos + len(seed)
            seqBm = seqBm[:endB]
            nuBm = nuBm[:, :endB]
            S = S[:endB]
            aliB = aliB[:, :endB]
            aliBm = aliBm[:, :endB]
        else:
            raise ValueError('ERROR: end not found')
    else:
        end = search_start + pos + len(seed)

    # If their lengths differ, align
    if (end - start) != len(seqBm):
        seqmt = seqm[start: end]
        from seqanpy import align_global
        ali = align_global(seqBm.tostring(), seqmt.tostring(), band=50)
        # FIXME: admit only gaps in our seq (gaps got stripped!)
        if '-' in ali[1]:
            raise ValueError('ERROR: gap in subtype B alignment (not implemented yet)')
        alim2 = np.fromstring(ali[2], 'S1')
        ind = alim2 != '-'
        seqBm = seqBm[ind]
        S = S[ind]
        nuBm = nuBm[:, ind]
        aliBm = aliBm[:, ind]

    return {'S': S, 'coord': (start, end),
            'seqB': seqBm.copy(),
            'nuB': nuBm.copy(),
            'aliBm': aliBm.copy(),
            'color': colors[protname]}


def get_entropy_pats(afts, VERBOSE=0):
    '''Calculate entropies'''
    n_patients = afts.shape[0]
    lseq = afts.shape[1]

    if VERBOSE >= 1:
        print 'Calculate entropy for each patient'
    S = -np.ones((n_patients, lseq), float)
    S_time = np.zeros((n_patients, lseq), object)
    times_covs = np.zeros(n_patients, object)
    for k, aft_pat in enumerate(afts):
        times_cov = aft_pat[0, 0] >= -0.1
        times_covs[k] = times_cov
        if times_cov.sum() < 2:
            S_time[k, :] = None
            continue
        
        for pos, aft_pos in enumerate(aft_pat):
            # FIXME: deal better with missing data (low coverage)
            Stmp = 0
            Stmp_time = np.zeros(times_cov.sum())
            for j, aft_nuc in enumerate(aft_pos):
                aft_nuc = aft_nuc[times_cov]
                Stmp -= ((aft_nuc + 1e-8) * np.log2(aft_nuc + 1e-8)).mean()
                Stmp_time -= ((aft_nuc + 1e-8) * np.log2(aft_nuc + 1e-8))
            S[k, pos] = Stmp
            S_time[k, pos] = Stmp_time

    return (S, S_time)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get shared allele trajectories',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the allele frequency trajectories')

    args = parser.parse_args()
    fragments = args.fragments
    VERBOSE = args.verbose
    plot = args.plot

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    plt.ioff()

    data_ali = {}

    for iif, fragment in enumerate(fragments):
        if VERBOSE >= 1:
            print fragment

        # Get our own coordinates, stripping gaps
        ali = AlignIO.read(root_patient_folder+'all/aft_shared_ali_'+fragment+'.fasta', 'fasta')
        if iif == 0:
            pnames = map(lambda x: x.name.split('_')[-1][1:], ali)
            patients = load_patients()
            patients = patients.loc[patients.index.isin(pnames)]

        alim = np.array(ali, 'S1')
        alim = alim[:, -((alim == '-').any(axis=0))]
        num = np.array([(alim == a).mean(axis=0) for a in alpha])
        seqm = alpha[num.argmax(axis=0)]

        # Low-coverage regions are bytecoded by -1
        from hivwholeseq.patients.filenames import root_patient_folder
        afts = np.load(root_patient_folder+'all/aft_shared_'+fragment+'.npy')
        n_patients = afts.shape[0]
        lseq = afts.shape[1]
        lalpha = afts.shape[2]

        if VERBOSE >= 1:
            print 'Get subtype B multiple sequence alignment'
        for protname in protnames[fragment]:
            if VERBOSE >= 1:
                print protname
            try:
                ddata = get_subtype_B_obs(protname, seqm, VERBOSE=VERBOSE)
                (start, end) = ddata['coord']
                ddata['alim'] = alim[:, start: end].copy()
                ddata['seqm'] = seqm[start: end].copy()
                # FIXME: we must use deepcopy here
                #ddata['afts'] = afts[:, start: end].copy()
                data_ali[protname] = ddata

            except ValueError:
                continue
 
        # Get entropies of patients
        (S, S_time) = get_entropy_pats(afts, VERBOSE=VERBOSE)
        for protname in protnames[fragment]:
            if protname not in data_ali:
                continue
            (start, end) = data_ali[protname]['coord']
            data_ali[protname]['S_pat'] = S[:, start: end].copy()
            data_ali[protname]['S_time'] = S_time[:, start: end].copy()       

        ## Look for parallel conservation
        #threshold = 0.005
        #n_poly = np.zeros((lseq, lalpha), int)  # Number of patients in which the allele is polymorphic at least once
        #n_cov = np.zeros(lseq, int)	                # Number of patients with sufficient coverage
        #for pos, aft_pos in enumerate(afts.swapaxes(0, 1).swapaxes(1, 2)): # Sort: pos, nuc, pats [, time] 
        #    if VERBOSE >= 3:
        #        print pos
        #    
        #    for j, aft_nuc in enumerate(aft_pos):
        #        n = 0
        #        nc = 0
        #        for aft_pat in aft_nuc:
        #            # FIXME: deal better with missing data (low coverage)
        #            times_cov = aft_pat >= -0.1
        #            if times_cov.sum() < 3:
        #                continue
        #            nc += 1
        #            aft_pat = aft_pat[times_cov]
        #            if ((aft_pat > threshold) & (aft_pat < 1 - threshold)).any():
        #                n += 1
        #        if j == 0:
        #            n_cov[pos] = nc
        #        n_poly[pos, j] = n
        #n_poly_max_raw = n_poly.max(axis=1)

        ## Exclude sites covered by few patients only
        #n_poly_max = n_poly_max_raw.copy()
        #n_poly_max[n_cov < 5] = -1

        ## Correlation coefficients between these things
        #corrs = defaultdict(lambda: defaultdict(list))
        #corrs_lowdiv = defaultdict(lambda: defaultdict(list))
        #corrs_highdiv = defaultdict(lambda: defaultdict(list))
        #Sthre = 0.1
        #for protname, ddata in data_ali.iteritems():
        #    (start, end) = ddata['coord']
        #    if VERBOSE >= 2:
        #        print protname
        #        print 'Spearman r between n_poly_max and S_B:',
        #        print spearmanr(n_poly_max_raw[start: end], ddata['S'])

        #        print 'Spearman r between <S>_pats and S_B:',
        #        print spearmanr(S[:, start: end].mean(axis=0), ddata['S'])

        #        print 'Spearman r between <S>_time and S_B, for each patient:'
        #        for i in xrange(n_patients):
        #            print pnames[i], spearmanr(S[i, start: end], ddata['S'])

        #        print ''
        #        print 'Spearman r between <S>_time and S_B, for each patient and time:'

        #    for i in xrange(n_patients):
        #        patient = Patient(patients.iloc[i])
        #        if S_time[i, 0] is None:
        #            continue
        #        for j in xrange(len(S_time[i, 0])):
        #            S_time_pat = np.array(map(itemgetter(j), S_time[i]))[start: end]
        #            S_subtypeB = ddata['S']
        #            (r, pvalue) = spearmanr(S_time_pat, S_subtypeB)
        #            corrs[protname][patient.name].append(r)

        #            ind_lowd = S_subtypeB < Sthre
        #            corrs_lowdiv[protname][patient.name].append(spearmanr(S_time_pat[ind_lowd], S_subtypeB[ind_lowd])[0])
        #            corrs_highdiv[protname][patient.name].append(spearmanr(S_time_pat[-ind_lowd], S_subtypeB[-ind_lowd])[0])
        #            if VERBOSE >= 2:
        #                print pnames[i], j, r, pvalue
        #        if VERBOSE >= 2:
        #            print ''

        #    if VERBOSE >= 2:
        #        print ''
        #    
        #if plot:
        #    if VERBOSE:
        #        print 'Plot along fragment'
        #    fig, ax = plt.subplots(figsize=(16, 5))
        #    ax.plot(np.arange(len(n_poly_max)), 1.0 * n_poly_max_raw / n_patients, lw=1.5, c='b', label='# pats')
        #    ax.plot(np.arange(len(n_poly_max)), 5.0 * S.mean(axis=0), lw=1.5, c='k', label='5 * <S>_pats')
        #    ax.plot(np.arange(len(n_poly_max)), 1.0 * n_cov / n_patients, lw=1.5, c='gray')
        #    for protname, ddata in data_ali.iteritems():
        #        (start, end) = ddata['coord']
        #        ax.plot(np.arange(start, end), ddata['S'], lw=1.5, c=ddata['color'], label='subtype B '+protname)
        #    ax.set_xlabel('Position [bp]')
        #    ax.set_ylabel('frac patients polymorphic')
        #    ax.set_title(fragment)
        #    ax.set_xlim(-1, len(n_poly_max) + 1)
        #    ax.set_ylim(-0.2, 1.9)
        #    ax.legend(loc=1)

        #    plt.tight_layout()


        #    if VERBOSE:
        #        print 'Plot correlation during infection'
        #    fig, axs = plt.subplots(1, len(corrs.keys()), figsize=(1 + len(corrs.keys()) * 7, 8))
        #    for i, (protname, corrs_prot) in enumerate(corrs.iteritems()):
        #        ax = axs[i]
        #        for j, (pname, corr) in enumerate(corrs_prot.iteritems()):
        #            color = cm.jet(1.0 * j / len(pnames))
        #            patient = Patient(patients.loc[pname])
        #            patient.discard_nonsequenced_samples
        #            times = (patient.times - patient.transmission_date) / np.timedelta64(1, 'D')
        #            times = times[times_covs[pnames.index(pname)]]
        #            ax.plot(times, corr, lw=2, label=pname, c=color)

        #        ax.set_xlabel('Time from transmission [days]')
        #        if i == 0:
        #            ax.set_ylabel('spearman r (S_sample, S_subtypeB)')
        #        ax.set_ylim(-0.02, 0.7)
        #        ax.set_title(protname)
        #        ax.legend(loc=4, fontsize=14)

        #    plt.tight_layout()


        #    if VERBOSE:
        #        print 'Plot correlation during infection stratified by conservation level'
        #    fig, axs = plt.subplots(1, len(corrs.keys()), figsize=(1 + len(corrs.keys()) * 7, 8))
        #    for i, protname in enumerate(corrs):
        #        ax = axs[i]
        #        for j, pname in enumerate(corrs_prot):
        #            corr_low = corrs_lowdiv[protname][pname]
        #            corr_high = corrs_highdiv[protname][pname]

        #            color = cm.jet(1.0 * j / len(pnames))
        #            patient = Patient(patients.loc[pname])
        #            patient.discard_nonsequenced_samples
        #            times = (patient.times - patient.transmission_date) / np.timedelta64(1, 'D')
        #            times = times[times_covs[pnames.index(pname)]]
        #            ax.plot(times, corr_low, lw=2, ls='-', label=pname+', S < '+str(Sthre), c=color)
        #            ax.plot(times, corr_high, lw=2, ls='--', label=pname+', S > '+str(Sthre), c=color)

        #        ax.set_xlabel('Time from transmission [days]')
        #        if i == 0:
        #            ax.set_ylabel('spearman r (S_sample, S_subtypeB)')
        #        ax.set_ylim(-0.02, 0.7)
        #        ax.set_title(protname)
        #        ax.legend(loc=4, fontsize=14)

        #    plt.tight_layout()

        #    
        #    if VERBOSE:
        #        print 'Plot entropy difference distribution with pat and time'
        #    for i in xrange(n_patients):
        #        S_pat = np.vstack(S_time[i])
        #        S_diff = {protname: (S_pat[d['coord'][0]: d['coord'][1]].T - d['S']).T
        #                  for protname, d in data_ali.iteritems()}

        #        S_diff = np.vstack(S_diff.itervalues())

        #        fig, ax = plt.subplots()
        #        pname = pnames[i]
        #        patient = Patient(patients.loc[pname])
        #        patient.discard_nonsequenced_samples
        #        times = (patient.times - patient.transmission_date) / np.timedelta64(1, 'D')
        #        times = times[times_covs[pnames.index(pname)]]
        #        for j in xrange(S_diff.shape[1]):
        #            Stmp = S_diff[:, j]
        #            Stmp.sort()
        #            x = Stmp
        #            y = np.linspace(1e-4, 1 - 1e-4, len(x))
        #            if use_logit:
        #                y = np.log10(y / (1 - y))
        #            ax.plot(x, y,
        #                    c=cm.jet(1.0 * j / S_diff.shape[1]),
        #                    label=str(times[j])+' days', lw=2)

        #        if use_logit:
        #            ax.set_ylim(-4.1, 4.1)
        #            trfun = lambda x: np.log10(x / (1 - x))
        #            tickloc = np.array([0.0001, 0.01, 0.5, 0.99, 0.9999])
        #            ax.set_yticks(trfun(tickloc))
        #            ax.set_yticklabels(map(str, tickloc))
        #            from matplotlib.ticker import FixedLocator
        #            ticklocminor = np.concatenate([[10**po * x for x in xrange(2 , 10)]
        #                                           for po in xrange(-4, -1)] + \
        #                                          [[0.1 * x for x in xrange(2 , 9)]] + \
        #                                          [[1 - 10**po * (10 - x) for x in xrange(2, 10)]
        #                                           for po in xrange(-2, -5, -1)])
        #            ax.yaxis.set_minor_locator(FixedLocator(trfun(ticklocminor)))
        #    
        #        ax.set_xlabel('S(pat, t) - S_subtypeB [bits]')
        #        ax.set_ylabel('Cum dist')
        #        ax.set_title(pnames[i])
        #        ax.legend(loc=2)
        #        ax.grid(True)

        #        plt.tight_layout()

    #if VERBOSE >= 2:
    #    print 'Most conserved sites in subtype B'
    #    print 'Prot  #   pos  |    S_B   pat <S_pat> | nB nC aft'
    #    print '---------------+----------------------+------------------------------------'
    #    for protname, d in data_ali.iteritems():
    #        ind_lowSB = np.argsort(d['S'])
    #        for jj, j in enumerate(ind_lowSB):
    #            if d['S'][j] > 0.01:
    #                continue
    #            for i, pname in enumerate(pnames):
    #                if d['S_pat'][i, j] > 0.05:
    #                    print '{:>3s}'.format(protname), '{:4d}'.format(jj), '{:4d}'.format(j),
    #                    print ' | ',
    #                    print '{:1.3f}'.format(d['S'][j]),
    #                    print '{:>6s}'.format(pname), '{:1.3f}'.format(d['S_pat'][i, j]),
    #                    print ' | ',
    #                    print d['seqB'][j], d['alim'][i, j],
    #                    print ' '.join(map('{:1.3f}'.format, d['afts'][i, j, alphal.index(d['alim'][i, j])]))

    #if VERBOSE >= 2:
    #    print 'Correlate initial allele and diversity'
    #    for protname, ddata in data_ali.iteritems():
    #        print protname
    #        for Slabel in ('SB', 'S_pat'):
    #            print Slabel
    #            print '   pat |  # B | # nB |  S mean B | S mean nB |  S median B |    S 25-75% B   | S median nB |    S 25-75% nB'
    #            print '-------+------+------+-----------+-----------+-------------+-----------------+-------------+----------------'
    #            for i, pname in enumerate(pnames):
    #                if Slabel == 'S_pat':
    #                    Stmp = ddata['S_pat'][i]
    #                else:
    #                    Stmp = ddata['S']
    #                seq_pat = ddata['alim'][i]
    #                seqB = ddata['seqB']
    #                ind_div = seq_pat != seqB
    #                print '{:6s}'.format(pname), '|',
    #                print '{:4d}'.format((-ind_div).sum()), '|',
    #                print '{:4d}'.format((+ind_div).sum()), '|',
    #                print '{:1.7f}'.format(Stmp[-ind_div].mean()), '|',
    #                print '{:1.7f}'.format(Stmp[+ind_div].mean()), '|',
    #                print '{:1.9f}'.format(np.median(Stmp[-ind_div])), '|',
    #                print ' '.join(map('{:1.1e}'.format, np.percentile(Stmp[-ind_div], (25, 75)))), '|',
    #                print '{:1.9f}'.format(np.median(Stmp[+ind_div])), '|',
    #                print ' '.join(map('{:1.1e}'.format, np.percentile(Stmp[+ind_div], (25, 75))))
    #            print ''

    #    print 'Stratify by subtype B entropy'
    #    n_strata = 7
    #    S_strata = np.zeros((2, len(pnames), n_strata), object) # NOTE: (B, nB) x pats x strata
    #    # Two columns: subtype B entropy and pat entropy
    #    for j in xrange(2):
    #        for i in xrange(len(pnames)):
    #            for k in xrange(n_strata):
    #                S_strata[j, i, k] = []
    #    for protname, ddata in data_ali.iteritems():
    #        print protname
    #        StmpB = ddata['S']
    #        ind_sort = np.argsort(StmpB)
    #        seqB = ddata['seqB']
    #        for i, pname in enumerate(pnames):
    #            seq_pat = ddata['alim'][i]
    #            ind_div = seq_pat != seqB
    #            Stmpp = ddata['S_pat'][i]

    #            for k in xrange(n_strata):
    #                # nB
    #                ind = np.intersect1d(ind_sort[k * len(ind_sort) // n_strata:
    #                                              (k + 1) * len(ind_sort) // n_strata],
    #                                     ind_div.nonzero()[0])
    #                Stmp = Stmpp[ind]
    #                Stmp = Stmp[-np.isnan(Stmp)]
    #                S_strata[1, i, k].extend(Stmp)

    #                # B
    #                ind = np.intersect1d(ind_sort[k * len(ind_sort) // n_strata:
    #                                              (k + 1) * len(ind_sort) // n_strata],
    #                                     (-ind_div).nonzero()[0])
    #                Stmp = Stmpp[ind]
    #                Stmp = Stmp[-np.isnan(Stmp)]
    #                S_strata[0, i, k].extend(Stmp)

    #    if plot:
    #        fig, ax = plt.subplots(1, 1)
    #        for k in xrange(n_strata):
    #            SnB = np.concatenate(S_strata[1, :, k])
    #            SnB.sort()
    #            SB = np.concatenate(S_strata[0, :, k])
    #            SB.sort()

    #            if k == 0:
    #                label = '1st quantile'
    #            elif k == 1:
    #                label = '2nd quantile'
    #            elif k == 2:
    #                label = '3rd quantile'
    #            else:
    #                label = str(k+1)+'th quantile'

    #            ax.plot(SB, np.linspace(0, 1, len(SB)), lw=2, c=cm.jet(1.0 * k / n_strata), ls='--', label=label+', B')
    #            if len(SnB) > 5:
    #                ax.plot(SnB, np.linspace(0, 1, len(SnB)), lw=2, c=cm.jet(1.0 * k / n_strata), ls='-', label=label+', nB')

    #        ax.set_xlabel('S pat [bits]')
    #        ax.set_ylabel('cum dist')
    #        ax.legend(loc=4, fontsize=10)
    #        ax.set_title('Entropy dist by stratum in subtype B')

    #        plt.tight_layout()

    #if plot:
    #    if VERBOSE:
    #        print 'Synonymous/nonsynonymous diversity'
    #    fig, axs = plt.subplots(1, len(data_ali), figsize=(2 + 5 * len(data_ali), 7))
    #    if len(data_ali) == 1:
    #        axs = [axs]
    #    from Bio.Seq import translate
    #    for ip, (protname, d) in enumerate(data_ali.iteritems()):
    #        seqBm = d['seqB']
    #        S_subtypeB = d['S']
    #        S_pat = d['S_pat']
    #        S_Bsns = {'syn': [], 'nonsyn': []}
    #        S_patsns = {'syn': [], 'nonsyn': []}
    #        for j in xrange(len(seqBm)):
    #            jcod = j // 3
    #            jp = j % 3
    #            cod_cons = seqBm[jcod: jcod + 3]
    #            cod_mut = cod_cons.copy()
    #            # Is there a possible syn mut at this site?
    #            n_deg = 0
    #            for a in alpha[:4]:
    #                cod_mut[jp] = a
    #                if translate(''.join(cod_cons)) == translate(''.join(cod_mut)):
    #                    n_deg += 1
    #            if n_deg >= 2:
    #                S_Bsns['syn'].append(S_subtypeB[j])
    #                S_patsns['syn'].append(S_pat[:, j])
    #            else:
    #                S_Bsns['nonsyn'].append(S_subtypeB[j])
    #                S_patsns['nonsyn'].append(S_pat[:, j])

    #        S_Bsns['syn'].sort()
    #        S_Bsns['nonsyn'].sort()
    #        S_patsns['syn'] = np.array(S_patsns['syn'])
    #        S_patsns['nonsyn'] = np.array(S_patsns['nonsyn'])
    #        S_patsns['syn'].sort(axis=1)
    #        S_patsns['nonsyn'].sort(axis=1)

    #        ax = axs[ip]
    #        ls = {'syn': '--', 'nonsyn': '-'}
    #        for key, S_tmp in S_Bsns.iteritems():
    #            ax.plot(S_tmp, np.linspace(0, 1, len(S_tmp)), label='B '+key, lw=2, ls=ls[key], c='grey')
    #            for i, pname in enumerate(pnames):
    #                S_tmp = S_patsns[key][i]
    #                ax.plot(5 * S_tmp, np.linspace(0, 1, len(S_tmp)), label=pname+' '+key, lw=2, ls=ls[key],
    #                        c=cm.jet(1.0 * i / len(pnames)))
    #        ax.set_xlabel('S_subtype [bits]')
    #        ax.set_title(protname)
    #        if ip == 0:
    #            ax.set_ylabel('cum dist')
    #            ax.legend(loc=4, fontsize=8)
    #        ax.set_xscale('log')
    #        ax.set_xlim(1e-3, 3)
    #        ax.grid(True)

    #    plt.tight_layout()

    if VERBOSE:
        print 'Correlation between single patient, few ones, and subtype'
        
        from itertools import combinations

        corrs = defaultdict(list)
        for protname, ddata in data_ali.iteritems():
            StmpB = ddata['S']
            for n in xrange(n_patients):
                corr = []
                for ind in combinations(xrange(n_patients), n + 1):
                    Stmps = ddata['S_pat'][list(ind)]
                    Stmps = Stmps.mean(axis=0)
                    r = spearmanr(StmpB, Stmps)[0]
                    corr.append(r)
                corrs[protname].append(corr)

        if plot:
            fig, axs = plt.subplots(1, len(data_ali), figsize=(2 + 5 * len(data_ali), 7))
            if len(data_ali) == 1:
                axs = [axs]
            iax = 0
            for protname in protnames_sort:
                if protname not in corrs:
                    continue
                
                corrprot = corrs[protname]
                ax = axs[iax]
                x = np.arange(n_patients) + 1
                y = corrprot
                ax.boxplot(y, positions=x)

                ax.set_ylim(0, 1)
                ax.set_xlabel('# patients')
                if iax == 0:
                    ax.set_ylabel('spearman r')
                else:
                    ax.set_yticklabels([])
                ax.set_title(protname)
                ax.grid(True)

                iax += 1

            fig.suptitle('Entropy correlation btw patients and subtype',
                         fontsize=18)

            plt.tight_layout(rect=(0, 0, 1, 0.93))


                    

    if plot:
        plt.ion()
        plt.show()

