# vim: fdm=marker
'''
author:     Fabio Zanini
date:       09/01/15
content:    Quantify purifying selection on different subtype entropy classes.
'''
# Modules
import os
from itertools import izip
from collections import defaultdict, Counter
import numpy as np
import pandas as pd

# FIXME: move to plots.py
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.paper_figures.filenames import get_figure_folder
from hivwholeseq.paper_figures.plots import plot_fitness_cost




# Globals
mu = 5e-6
pnames = ['20097', '15376', '15823', '15241', '9669', '15319']
regions = ['p17', 'p24', 'nef', 'PR', 'RT', 'IN', 'vif']
fun = lambda x, l, u: l * (1 - np.exp(- u/l * x))



# Functions
def fit_fitness_cost(x, y, VERBOSE=0, mu=5e-6):
    '''Fit saturation curve for fitness costs
    
    NOTE: as can be seen from the plots below, the fit for l is not very
    sensitive to mu.
    '''
    fun_min = lambda x, l, u: ((y - fun(x, l, u))**2).sum()

    # Minimization poor man's style: we compute the function everywhere
    # NOTE: I am not retarded. This solution is chosen because it's still cheap
    # and the right flank of the objective function in l is much steeper than
    # the left flank, causing some troubles
    ls = np.logspace(-5, 0)
    z = np.array([fun_min(x, l, mu) for l in ls])
    l = ls[z.argmin()]

    # FIXME: if we fix mu, this is useless, but it's ok for now
    params = (l, mu)

    return params


def plot_function_minimization(x, y, params):
    '''Investigate inconsistencies in fits'''
    fun_min = lambda p: ((y - fun(x, p[0], p[1]))**2).sum()

    p1 = np.logspace(np.log10(params[0]) - 3, np.log10(params[0]) + 3, 10)
    p2 = np.logspace(np.log10(params[1]) - 3, np.log10(params[1]) + 3, 10)

    p1G = np.tile(p1, (len(p2), 1))
    p2G = np.tile(p2, (len(p1), 1)).T
    pG = np.dstack([p1G, p2G])
    z = np.log(np.array([[fun_min(ppp) for ppp in pp] for pp in pG]))

    fig, ax = plt.subplots()
    ax.imshow(z, interpolation='nearest')

    ax.set_xlabel('Log l')
    ax.set_ylabel('Log u')

    plt.ion()
    plt.show()


def plot_function_minimization_1d(x, y, l, us=[1.2e-6], title=''):
    '''Investigate inconsistencies in fits'''
    fun_min = lambda l, u: ((y - fun(x, l, u))**2).sum()

    p1 = np.logspace(np.log10(l) - 3, np.log10(l) + 3, 100)
    zs = np.log(np.array([[fun_min(pp, u) for pp in p1] for u in us]))

    fig, ax = plt.subplots()

    from itertools import izip
    for i, (z, u) in enumerate(izip(zs, us)):
        ax.plot(p1, z, lw=2, color=cm.jet(1.0 * i / len(us)),
                label='mu = {:1.1e}'.format(u))

    if title:
        ax.set_title(title)
    ax.set_xlabel('Saturation frequency ($\mu / s$)')
    ax.set_xscale('log')
    ax.grid(True)
    ax.legend(loc='upper left')

    plt.ion()
    plt.show()


def fit_purifying(pnames, regions, VERBOSE=0):
    '''Collect fits for all regions'''
    from hivwholeseq.miseq import alpha, alphal
    from hivwholeseq.patients.patients import load_patients, Patient
    from hivwholeseq.analysis.explore_entropy_patsubtype import (
        get_subtype_reference_alignment, get_ali_entropy)
    from hivwholeseq.cross_sectional.get_subtype_entropy import (
        get_subtype_reference_alignment_entropy)
    from hivwholeseq.cross_sectional.get_subtype_consensus import (
        get_subtype_reference_alignment_consensus)
    from hivwholeseq.patients.one_site_statistics import get_codons_n_polymorphic
    from hivwholeseq.analysis.mutation_rate.explore_divergence_synonymous import translate_masked

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

            # FIXME: deal better with depth (this should be already there?)
            aft[aft < 2e-3] = 0

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

    # Bin by subtype entropy, using quantiles
    #bins_S = np.array([0, 0.03, 0.06, 0.1, 0.25, 0.7, 3])
    bins_S = np.array(data['Ssub'].quantile(q=np.linspace(0, 1, 6)))
    binsc_S = 0.5 * (bins_S[1:] + bins_S[:-1])
    data['Sbin'] = 0
    for b in bins_S[1:]:
        data.loc[data.loc[:, 'Ssub'] >= b, 'Sbin'] += 1


    # Fit exponential saturation
    mu = 5e-6
    fits = []
    dataf = (data
             .loc[data.loc[:, 'Ssub'] < bins_S[-1]]
             .loc[:, ['region', 'Sbin', 'time', 'af']]
             .groupby(['region', 'Sbin']))
    for (region, iSbin), datum in dataf:
        x = np.array(datum['time'])
        y = np.array(datum['af'])

        ind = -(np.isnan(x) | np.isnan(y))
        x = x[ind]
        y = y[ind]

        try:
            (l, u) = fit_fitness_cost(x, y, mu=mu)
            if VERBOSE >= 3:
                plot_function_minimization_1d(x, y, l, us=[1e-6, 2e-6, 5e-6, 1e-5],
                                              title=region+', iSbin = '+str(iSbin))

        except RuntimeError:
            continue

        fits.append((region, iSbin, l, u))

    fits = pd.DataFrame(data=fits,
                        columns=['region', 'iSbin', 'l', 'u'])
    fits['S'] = binsc_S[fits['iSbin']]
    fits['Smin'] = bins_S[fits['iSbin']]
    fits['Smax'] = bins_S[fits['iSbin'] + 1]
    
    # Estimate fitness cost
    fits['s'] = mu / fits['l']

    return fits


def plot_single_fit(region, fitsreg, title='', VERBOSE=0):
    '''Plot the fits for purifying selection'''

    fig, ax = plt.subplots(figsize=(6, 6))
    xfit = np.logspace(0, 3.5, 1000)

    for _, fit in fitsreg.iterrows():
        iSbin = fit['iSbin']
        Smin = fit['Smin']
        Smax = fit['Smax']
        l = fit['l']
        u = fit['u']
        yfit = fun(xfit, l, u)
        label = ('S e ['+'{:2.2f}'.format(Smin)+', '+'{:2.2f}'.format(Smax)+']'+
                 ', s = '+'{:.1G}'.format(mu / l))

        
        color = cm.jet(1.0 * iSbin / len(fitsreg))

        ax.plot(xfit, yfit, color=color, label=label, lw=2)

    ax.set_xlabel('Time [days from infection]')
    ax.set_ylabel('Allele frequency')
    ax.legend(loc='lower right', title='Entropy class', fontsize=10)
    ax.text(0.05, 0.9,
            ('$f(t) \, = \, \mu / s \, [1 - e^{-st}]$'),
            fontsize=20,
            horizontalalignment='left',
            verticalalignment='center',
            transform=ax.transAxes)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True)



# Script
if __name__ == '__main__':

    VERBOSE = 2
    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'
    mkdirs(fn_data)
    fn_data = fn_data + 'fitness_cost.pickle'

    if not os.path.isfile(fn_data):
        fits = fit_purifying(pnames, regions, VERBOSE=VERBOSE)

        if VERBOSE >= 1:
            print 'Save to file'
        fits.to_pickle(fn_data)

    else:
        fits = pd.read_pickle(fn_data)

    ## Plot the fits, one fig per region
    #for (region, fitsreg) in fits.groupby('region'):
    #    plot_single_fit(region, fitsreg)

    #filename = foldername+'fitness_cost_by_region'
    #for ext in ['png', 'pdf', 'svg']:
    #    plot_fitness_cost(fits,
    #                      VERBOSE=VERBOSE,
    #                      savefig=filename+'.'+ext)

    plot_fitness_cost(fits, VERBOSE=VERBOSE)

    # TODO: add bootstrap

