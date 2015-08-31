# vim: fdm=marker
'''
author:     Fabio Zanini
date:       09/01/15
content:    Quantify purifying selection on different subtype entropy classes,
            and, at the same time, the mutation rate.
'''
# Modules
import os
import sys
import argparse
from itertools import izip
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha, alphal
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.utils.sequence import translate_with_gaps
import hivwholeseq.utils.plot
from hivwholeseq.analysis.explore_entropy_patsubtype import (
    get_subtype_reference_alignment, get_ali_entropy)
from hivwholeseq.cross_sectional.get_subtype_entropy import (
    get_subtype_reference_alignment_entropy)
from hivwholeseq.cross_sectional.get_subtype_consensus import (
    get_subtype_reference_alignment_consensus)
from hivwholeseq.patients.one_site_statistics import get_codons_n_polymorphic
from hivwholeseq.analysis.mutation_rate.explore_divergence_synonymous import translate_masked
from hivwholeseq.analysis.purifying_selection.filenames import get_fitness_cost_entropy_filename

from hivwholeseq.analysis.purifying_selection.quantify_selection_genomewide import (
    plot_fits as plot_fitness_cost)
from hivwholeseq.analysis.mutation_rate.comparison_Abram import plot_comparison
from hivwholeseq.analysis.mutation_rate.mutation_rate import plot_mu_matrix



# Globals
pnames = ['p1', 'p3', 'p5', 'p8', 'p9', 'p11']
regions = ['p17', 'p24', 'nef', 'PR', 'RT', 'IN', 'vif']
fun = lambda x, l, u: l * (1 - np.exp(- u/l * x))



# Functions
def add_Sbins(data, bins=8, VERBOSE=0):
    '''Add entropy bins to the data'''
    if np.isscalar(bins):
        bins = np.array(data['Ssub'].quantile(q=np.linspace(0, 1, bins + 1)))
    
    binsc = 0.5 * (bins[1:] + bins[:-1])

    data['Sbin'] = -1
    for b in bins:
        data.loc[data.loc[:, 'Ssub'] >= b, 'Sbin'] += 1
    data['Sbin'] = data['Sbin'].clip(0, len(binsc) - 1)

    return bins, binsc


def collect_data_joint(pnames, regions, VERBOSE=0, plot=False):
    '''Collect data for joint estimate of fitness cost and mutation rates'''

    data = []

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    for region in regions:
        if VERBOSE >= 1:
            print region

        if VERBOSE >= 2:
            print 'Get subtype consensus (for checks only)'
        conssub = get_subtype_reference_alignment_consensus(region, VERBOSE=VERBOSE)

        if VERBOSE >= 2:
            print 'Get subtype entropy'
        Ssub = get_subtype_reference_alignment_entropy(region, VERBOSE=VERBOSE)

        for ipat, (pname, patient) in enumerate(patients.iterrows()):
            pcode = patient.code
            if VERBOSE >= 2:
                print pname, pcode

            patient = Patient(patient)
            aft, ind = patient.get_allele_frequency_trajectories(region,
                                                                 cov_min=1000,
                                                                 depth_min=300,
                                                                 VERBOSE=VERBOSE)
            if len(ind) == 0:
                if VERBOSE >= 2:
                    print 'No time points: skip'
                continue

            times = patient.times[ind]

            if VERBOSE >= 2:
                print 'Get coordinate map'
            coomap = patient.get_map_coordinates_reference(region, refname=('HXB2', region))

            icons = patient.get_initial_consensus_noinsertions(aft, VERBOSE=VERBOSE,
                                                               return_ind=True)
            consm = alpha[icons]
            protm = translate_masked(consm)
            
            # Premature stops in the initial consensus???
            if '*' in protm:
                # Trim the stop codon if still there (some proteins are also end of translation)
                if protm[-1] == '*':
                    if VERBOSE >= 2:
                        print 'Ends with a stop, trim it'
                    icons = icons[:-3]
                    consm = consm[:-3]
                    protm = protm[:-1]
                    aft = aft[:, :, :-3]
                    coomap = coomap[coomap[:, 1] < len(consm)]

                else:
                    continue

            # Get the map as a dictionary from patient to subtype
            coomapd = {'pat_to_subtype': dict(coomap[:, ::-1]),
                       'subtype_to_pat': dict(coomap)}

            # Get only codons with at most one polymorphic site, to avoid obvious epistasis
            ind_poly, _ = get_codons_n_polymorphic(aft, icons, n=[0, 1], VERBOSE=VERBOSE)
            ind_poly_dna = [i * 3 + j for i in ind_poly for j in xrange(3)]

            for posdna in ind_poly_dna:
                # Get the entropy
                if posdna not in coomapd['pat_to_subtype']:
                    continue
                pos_sub = coomapd['pat_to_subtype'][posdna]
                if pos_sub >= len(Ssub):
                    continue
                Ssubpos = Ssub[pos_sub]

                # Get allele frequency trajectory
                aftpos = aft[:, :, posdna].T

                # Get only non-masked time points
                indpost = -aftpos[0].mask
                if indpost.sum() == 0:
                    continue
                timespos = times[indpost]
                aftpos = aftpos[:, indpost]

                anc = consm[posdna]
                ianc = icons[posdna]

                # Skip if the site is already polymorphic at the start
                if aftpos[ianc, 0] < 0.95:
                    continue

                # Skip if the site has sweeps (we are looking at purifying selection only)
                # Obviously, it is hard to distinguish between sweeps and unconstrained positions
                # but that's not crucial because sweeps happen at unconstrained positions anyway
                if (aftpos[ianc] < 0.6).any():
                    continue

                for inuc, af in enumerate(aftpos[:4]):
                    nuc = alpha[inuc]
                    if nuc == anc:
                        continue

                    mut = anc+'->'+nuc

                    # Define transition/transversion
                    if frozenset(nuc+anc) in (frozenset('CT'), frozenset('AG')):
                        trclass = 'ts'
                    else:
                        trclass = 'tv'

                    # Get the whole trajectory for plots against time
                    for af, time in izip(aftpos[inuc], timespos):
                        data.append((region, pcode,
                                     posdna, pos_sub,
                                     anc, nuc, mut,
                                     trclass,
                                     Ssubpos,
                                     time, af))

    data = pd.DataFrame(data=data,
                        columns=['region', 'pcode',
                                 'posdna', 'possub',
                                 'anc', 'der', 'mut',
                                 'tr',
                                 'Ssub',
                                 'time', 'af'])

    return data


def fit_joint_model(data, method='group', VERBOSE=0, plot=False):
    '''Fit fitness cost and mutation rates at once'''
    from scipy.optimize import minimize

    def get_initial_parameters(muts):
        '''Get initial parameter vector for the fit'''
    
        # Get initial estimate for mu from Abram2010, resorted, PER DAY
        from hivwholeseq.analysis.mutation_rate.comparison_Abram import get_mu_Abram2010
        mu0 = get_mu_Abram2010().loc[muts] / 2.0
    
        # Get initial estimate for s from our fixed-mu estimate
        from hivwholeseq.analysis.purifying_selection.filenames import get_fitness_cost_entropy_filename
        sdata = pd.read_pickle(get_fitness_cost_entropy_filename('all_Abram2010'))
        s = sdata['s']
    
        p0 = np.concatenate(map(np.array, [mu0, s]))
        return {'p0': p0, 'mu0': mu0, 's0': s, 'sdata': sdata}

    # NOTE: muts just sets the order of mutations in the vector
    muts = np.unique(data['mut']).tolist()
    pdict = get_initial_parameters(muts)
    p0 = pdict['p0']
    s0 = pdict['s0']
    mu0 = pdict['mu0']
    sdata = pdict['sdata']

    # Bin by mutation type and subtype entropy, taking bins from fitness cost estimate
    #FIXME: we could also make time bins (less data, more smooth)
    data['mbin'] = map(muts.index, data['mut'])
    bins_S = np.insert(np.array(sdata['Smax']), 0, 0)
    bins_S, binsc_S = add_Sbins(data, bins=bins_S, VERBOSE=VERBOSE)

    if method == 'single':
        dataf = (data
                 .loc[:, ['mbin', 'Sbin', 'time', 'af']]
                 .groupby(['mbin', 'Sbin']))
    else:
        dataf = (data
                 .loc[:, ['mbin', 'Sbin', 'time', 'af']]
                 .groupby(['mbin', 'Sbin', 'time'])
                 .mean())
        dataf['time'] = dataf.index.get_level_values('time')
        dataf['Sbin'] = dataf.index.get_level_values('Sbin')
        dataf['mbin'] = dataf.index.get_level_values('mbin')
        dataf = dataf.groupby(['mbin', 'Sbin'])

    def fun_min(p):
        '''p is the parameters vector: #I mut rates + #J selection coefficients'''
        fun = lambda x, u, s: u / s * (1 - np.exp(- s * x))

        mu = p[:len(muts)]
        s = p[len(muts):]

        #res = 0
        #for (imb, isb), datum in dataf:
        #    res += ((datum['af'] - fun(datum['time'], mu[imb], s[isb]))**2).sum()

        res = sum(((datum['af'] - fun(datum['time'], mu[imb], s[isb]))**2).sum()
                  for (imb, isb), datum in dataf)

        return res
        

    def get_funmin_neighbourhood(fun_min, p0, nfev=30, plot=False, title=''):
        '''Calculate a few function evaluations for testing'''
        dynrangeexp = 0.5
        factorl = np.sort(np.concatenate([np.logspace(-dynrangeexp, dynrangeexp, nfev - 1),
                                         [1]]))

        datap = []
        for i in xrange(len(p0)):
            if VERBOSE >= 2:
                print 'Parameter n.'+str(i+1)
            x = []
            y = []

            for factor in factorl:
                p = p0.copy()
                p[i] *= factor

                x.append(factor)
                y.append(fun_min(p))

            if i < len(muts):
                label = muts[i]
            else:
                label = 'S ~ {:.1G}'.format(binsc_S[i - len(muts)])

            color = cm.jet(1.0 * i / len(p0))

            datap.append({'x': x, 'y': y, 'label': label, 'color': color})

        # Plot
        if plot:
            fig, ax = plt.subplots()
            for datum in datap:
                x = datum['x']
                y = datum['y']
                label = datum['label']
                color = datum['color']

                ax.plot(x, y, lw=2, label=label, color=color)

            ax.set_xlabel('Factor change of parameter')
            ax.set_xscale('log')
            ax.set_ylabel('Fmin')
            ax.set_yscale('log')
            ax.set_xlim(10**(-dynrangeexp), 10**dynrangeexp)

            ax.grid(True)
            ax.legend(loc='upper center', ncol=3, fontsize=8)

            plt.tight_layout()

        return datap

    if plot:
        if VERBOSE >= 1:
            print 'Plot objective function in a neighbourhood of initial parameters'

        datap0 = get_funmin_neighbourhood(fun_min, p0, plot=plot)

    if VERBOSE >= 1:
        print 'Minimize joint objective function'
    
    method = 'Powell' # minimize along each parameter at every iteration
    res = minimize(fun_min, p0, method=method,
                   options={'disp': True, 'maxiter': 10})
    print res
    pmin = res.x
    
    if plot:
        if VERBOSE >= 1:
            print 'Plot objective function in a neighbourhood of the solution'

        datapmin = get_funmin_neighbourhood(fun_min, pmin, plot=plot)

        plt.ion()
        plt.show()

    # Save results of minimization
    mu_min = pd.Series(pmin[:len(muts)], index=muts)
    s_min = pd.DataFrame(pmin[len(muts):], columns=['s'])
    s_min['Smin'] = sdata['Smin']
    s_min['Smax'] = sdata['Smax']
    s_min['S'] = sdata['S']
    s_min['iSbin'] = s_min.index

    return mu_min, s_min


def plot_mutation_rate(mu, VERBOSE=0):
    '''Plot mutation rate after joint minimization'''

    # Plot matrix
    muM = np.ma.masked_all((4, 4))
    for ia1, a1 in enumerate(alpha[:4]):
        for ia2, a2 in enumerate(alpha[:4]):
            if a1+'->'+a2 in mu.index:
                muM[ia1, ia2] = mu.loc[a1+'->'+a2]

    fig, ax = plt.subplots()
    h = ax.imshow(np.log10(muM.T + 1e-10),
                  interpolation='nearest',
                  vmin=-9, vmax=-4)

    ax.set_xticks(np.arange(4))
    ax.set_yticks(np.arange(4))
    ax.set_xticklabels(alpha[:4])
    ax.set_yticklabels(alpha[:4])

    ax.set_xlabel('From:')
    ax.set_ylabel('To:')

    cb = fig.colorbar(h)
    cb.set_ticks(np.arange(-9, -3))
    cb.set_ticklabels(['$10^{'+str(x)+'}$' for x in xrange(-9, -3)])
    cb.set_label('changes / position / day', rotation=270, labelpad=30)

    plt.tight_layout(rect=(-0.1, 0, 1, 1))

    # Plot comparison with Abram2010
    from hivwholeseq.analysis.mutation_rate.comparison_Abram import add_Abram2010
    mus = add_Abram2010(mu, VERBOSE=VERBOSE)

    fig, ax = plt.subplots()
    x = np.array(mus['Abram2010'])
    y = np.array(mus['new'])

    r = np.corrcoef(np.log10(x), np.log10(y))[0, 1]

    ax.scatter(x, y, s=40, color='k', label='{:2.0%}'.format(r))
    xl = np.logspace(-8, -4, 100)
    ax.plot(xl, xl, lw=2, c='grey')
    ax.set_xlabel('Abram et al 2010 [changes / site / day]')
    ax.set_ylabel('Joint estimate [changes / site / day]')
    ax.set_xlim(1e-8, 1e-4)
    ax.set_ylim(1e-8, 1e-4)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='upper left', title='Correlation\ncoefficient:', fontsize=14)
    ax.grid(True)

    plt.tight_layout()



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(
        description='Infer mutation rates AND fitness costs',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+', default=pnames,
                        help='Patient to analyze')
    parser.add_argument('--regions', nargs='+', default=regions,
                        help='Regions to analyze (e.g. F1 p17)')
    parser.add_argument('--verbose', type=int, default=2,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot results')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    plot = args.plot

    data = collect_data_joint(pnames, regions, VERBOSE=VERBOSE, plot=plot)

    mu_min, s_min = fit_joint_model(data, VERBOSE=VERBOSE, plot=plot)    

    if plot:
        data['mu'] = np.array(mu_min.loc[data['mut']])
        plot_fitness_cost(s_min, VERBOSE=VERBOSE, data=data)

        plot_mu_matrix(mu_min, time_unit='days')
        plot_comparison(mu_min)

        plt.ion()
        plt.show()

    if VERBOSE >= 1:
        print 'Save to file'
    from hivwholeseq.analysis.filenames import analysis_data_folder
    fn_out_mu = analysis_data_folder+'mu_joint.pickle'
    mu_min.to_pickle(fn_out_mu)

    from hivwholeseq.analysis.purifying_selection.filenames import get_fitness_cost_entropy_filename
    fn_out_s = get_fitness_cost_entropy_filename('all_joint')
    s_min.to_pickle(fn_out_s)

