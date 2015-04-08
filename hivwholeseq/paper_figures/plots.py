# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/12/14
content:    Module collecting all plot functions for unity of themes and similia.
'''
# Globals
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm

from hivwholeseq.utils.generic import mkdirs


# Functions
def plot_cuts_quality_along_reads(data, title='',
                                  VERBOSE=0, savefig=False):
    '''Plot some cuts of the quality along the read'''
    plt.ioff()

    from itertools import izip

    if VERBOSE:
        print 'Plot quality along read'

    fig, axs = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    axs[0].set_ylabel('Percentage of bases above quality x', fontsize=14)

    sns.set_style('whitegrid')
    colors = sns.dark_palette("blue", input="xkcd", n_colors=len(data[0]))

    for i, (ax, data_read) in enumerate(izip(axs, data)):
        for j, datum in enumerate(data_read):
            x = datum['x']
            y = datum['y']
            qthresh = datum['threshold']
            ax.plot(x, y,
                    color=colors[j],
                    lw=2,
                    label='Q = '+str(qthresh))
        ax.set_xlabel('Position [bp]', fontsize=14)
        ax.set_title('Read'+str(i+1), fontsize=16)
        ax.set_ylim(-1, 101)
        ax.set_xlim(-1, len(x) + 1)
        ax.legend(loc='lower center')
        ax.grid(True)

    if title:
        fig.suptitle(title, fontsize=20)

    if savefig:
        fig_filename = savefig
        fig_folder = os.path.dirname(fig_filename)

        mkdirs(fig_folder)
        fig.savefig(fig_filename)
        plt.close(fig)

    else:
        plt.tight_layout()
        plt.ion()
        plt.show()


def plot_phred_errors(data, label='', title='', VERBOSE=0, savefig=False):
    '''Plot the error counts'''
    plt.ioff()

    if VERBOSE:
        print 'Plot phred errors'

    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    sns.set_style('whitegrid')

    for i, datum in enumerate(data):
        counts = datum['counts']
        label = datum['label']

        y = counts.sum(axis=1).sum(axis=0)
        y = 1.0 * y[:, 1] / y[:, 0]
        x = np.arange(len(y))
        ax.plot(x, y, lw=2, alpha=0.8)
        ax.scatter(x, y, lw=2, label=label)
    
    ax.set_xlabel('Phred score')
    ax.set_ylabel('Error rate')
    ax.set_yscale('log')
    ax.set_ylim(1e-5, 1)
    ax.set_xlim(0, 45)
    ax.grid(True)

    if title:
        ax.set_title(title)

    # Reference line
    ax.plot(np.arange(len(y)), 10**(-0.1 * np.arange(len(y))), lw=1,
            color='k', ls='--')
    
    if savefig:
        fig_filename = savefig
        fig_folder = os.path.dirname(fig_filename)

        mkdirs(fig_folder)
        fig.savefig(fig_filename)
        plt.close(fig)

    else:
        plt.tight_layout()
        plt.ion()
        plt.show()
 

def plot_minor_allele_reference(data, title='', VERBOSE=0, savefig=False):
    '''Plot minor allele in the reference sample'''
    plt.ioff()

    if VERBOSE:
        print 'Plot minor alleles of control sample'

    fig, axs = plt.subplots(1, 2, figsize=(8, 5), sharey=True,
                            gridspec_kw={'width_ratios': [3, 1]})
    sns.set_style('whitegrid')

    for datum in data:
        y = datum['freq_minor']
        x = np.arange(len(y))
        axs[0].plot(x, y, lw=1.5, alpha=0.8)
        axs[0].scatter(x, y, lw=1.5)

        h = np.histogram(y, bins=np.logspace(-6, 0, 40))
        axs[1].barh(h[1][:-1], h[0], (h[1][1:] - h[1][:-1]))

    axs[0].set_xlabel('Position [bp]')
    axs[0].set_ylabel('Minor allele frequency (errors)')
    axs[0].set_yscale('log')
    axs[0].set_ylim(1e-6, 1)
    axs[0].set_xlim(-5, len(x) + 5)
    axs[0].grid(True)

    axs[1].set_xlabel('Number of alleles')
    axs[1].grid(True)
    axs[1].set_yscale('log')
    axs[1].set_xlim(0, 1.15 * h[0].max())

    if title:
        fig.suptitle(title)

    if savefig:
        fig_filename = savefig
        fig_folder = os.path.dirname(fig_filename)

        mkdirs(fig_folder)
        fig.savefig(fig_filename)
        plt.close(fig)

    else:
        plt.tight_layout()
        plt.ion()
        plt.show()


def plot_template_numbers(data, VERBOSE=0, title='', savefig=False):
    '''Plot template numbers'''
    plt.ioff()

    if VERBOSE:
        print 'Plot template numbers'

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    sns.set_style('whitegrid')

    ax = axs[0]
    agemax = np.concatenate([data_pat['age'] for data_pat in data]).max()
    nmax = 1e5

    # Plot reference (PCR eta = 1)
    ax.plot([1, nmax], [1, nmax], lw=2)

    colors = sns.dark_palette("cerulean", n_colors=100, reverse=True, input='xkcd')
    for i, data_pat in enumerate(data):
        ci = np.array(len(colors) * data_pat['age'] / agemax, int)
        ci[ci >= len(colors)] = len(colors) - 1
        color = colors[ci]
        ax.scatter(data_pat['n_approx'],
                   np.array(data_pat['n_dil']) * (1.05 * i / len(data)),
                   s=30,
                   color=color)
    
    ax.set_xlabel('# templates (from viral load)')
    ax.set_ylabel('# templates (from dilutions)')
    ax.set_ylim((1, nmax))
    ax.set_xlim((1, nmax))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True)

    ax = axs[1]
    x = np.concatenate([data_pat['n_approx'] for data_pat in data])
    y = np.concatenate([data_pat['n_dil'] for data_pat in data])
    ind = -(np.isnan(x) | np.isnan(y) | (x == 0) | (y == 0))
    x = x[ind]
    y = y[ind]
    age = np.concatenate([data_pat['age'] for data_pat in data])[ind]

    h = np.histogram(y / x, bins=np.logspace(-3, 0, 10), density=False)
    ax.bar(h[1][:-1], h[0], width=np.diff(h[1]))
    ax.set_xlabel('PCR efficiency')
    ax.set_ylabel('# of samples')
    ax.grid(True)

    ax.set_xscale('log')
    
    if title:
        fig.suptitle(title)

    if savefig:
        fig_filename = savefig
        fig_folder = os.path.dirname(fig_filename)

        mkdirs(fig_folder)
        fig.savefig(fig_filename)
        plt.close(fig)

    else:
        plt.tight_layout()
        plt.ion()
        plt.show()


def plot_insert_size_distribution(data, title='', VERBOSE=0, savefig=False):
    '''Plot histogram of insert sizes'''
    plt.ioff()

    if VERBOSE:
        print 'Plot template numbers'
    
    fig, ax = plt.subplots(1, 1)
    sns.set_style('whitegrid')

    colors = sns.xkcd_palette(['rusty red', 'slate blue'])
    n_data = len(data) + 0.1
    for i, datum in enumerate(data):
        w = np.diff(datum['bins'])
        x_offset = 1.0 * i / n_data * w
        w = 1.0 * w / n_data
        ax.bar(datum['bins'][:-1] + x_offset, datum['counts'], w,
               color=colors[i], label=datum['label'])

    ax.set_xlabel('Insert size')
    ax.set_ylabel('Number of inserts')
    ax.set_ylim(0, ax.get_ylim()[1] * 1.15)
    ax.legend(loc=1)

    if title:
        ax.set_title(title)

    if savefig:
        fig_filename = savefig
        fig_folder = os.path.dirname(fig_filename)

        mkdirs(fig_folder)
        fig.savefig(fig_filename)
        plt.close(fig)

    else:
        plt.tight_layout()
        plt.ion()
        plt.show()


def plot_rna_recombination(data, title='', VERBOSE=0, savefig=False):
    '''Plot recombination in the RNA mixes'''
    plt.ioff()

    if VERBOSE >= 1:
        print 'Plot RNA mix amplification bias and conditional switch probability'

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    sns.set_style('whitegrid')

    for datum in data:
        switchesn = datum['switchesn']
        counts = datum['counts']
        refnames = datum['refnames']
        al_polyd = datum['aldict']
        L = np.max(counts.keys())

        # Amplification bias
        ax = axs[0]
        ax.set_xlabel('Position [bp]')
        ax.set_ylabel('Frequency of '+refnames[0])
        data = [(pos, 1.0 * co[al_polyd[pos][1]] / sum(co.itervalues()))
                for (pos, co) in counts.iteritems() if sum(co.itervalues()) > 100]
        (x, y) = np.array(data).T
        ind = x.argsort()
        x = x[ind]
        y = y[ind]
        y[y > 1-1e-4] = 1-1e-4
        y[y < 1e-4] = 1e-4
        ax.plot(x, y, lw=2)
        ax.set_ylim(0.9e-4, 1 - 0.9e-4)
        ax.set_yscale('logit')
        ax.grid(True)

        # Recombination rate
        ax = axs[1]
        ax.set_xlabel('Distance [bp]')
        ax.set_ylabel('Conditional switch probability')
    
        colors = sns.dark_palette("cerulean", n_colors=100, reverse=True, input='xkcd')
        for (pos1, pos2), freq in switchesn.iteritems():
            posm = 0.5 * (pos1 + pos2)
            ci = np.array(len(colors) * posm / L, int)
            ci[ci >= len(colors)] = len(colors) - 1
            color = colors[ci]
            ax.scatter(pos2 - pos1, freq + 1e-6, s=30, color=color)

        ax.set_ylim(-0.05, 1.1)
        ax.grid(True)

    if title:
        fig.suptitle(title)

    if savefig:
        fig_filename = savefig
        fig_folder = os.path.dirname(fig_filename)

        mkdirs(fig_folder)
        fig.savefig(fig_filename)
        plt.close(fig)

    else:
        plt.tight_layout()
        plt.ion()
        plt.show()



def plot_allele_frequency_overlap(data, title='', VERBOSE=0, use_logit=True, savefig=False):
    '''Plot allele frequency in the overlap regions'''

    fig, ax = plt.subplots()
    sns.set_style('whitegrid')

    colors = sns.color_palette('Set1', 5)
    legend = set()
    for ida, datum in enumerate(data):
        afjoint = datum['af']
        color = colors[datum['io']]
        if datum['overlap'] not in legend:
            label = datum['overlap']
            legend.add(datum['overlap'])
        else:
            label = None
        ax.scatter(afjoint[0].ravel(), afjoint[1].ravel(),
                   s=50,
                   color=color,
                   alpha=0.7,
                   label=label,
                   edgecolor='none')
    
    xmin = 1e-4
    ax.plot([xmin, 1 - xmin], [xmin, 1 - xmin], lw=2, color='k', alpha=0.5)

    ax.set_xlabel('allele frequency leading fragment')
    ax.set_ylabel('allele frequency trailing fragment')
    ax.grid(True)

    if use_logit:
        ax.set_xscale('logit')
        ax.set_yscale('logit')
        ax.set_xlim(xmin, 1 - xmin)
        ax.set_ylim(xmin, 1 - xmin)

    # Plot stddev of a certain number of molecules in Poisson sampling
    n = 300
    x = np.linspace(-4, 0, 1000)
    x = 1.0 / (1 + 10**(-x))
    y = x - np.sqrt(x / n)
    ax.plot(np.concatenate([x, 1 - y[::-1]]), np.concatenate([y, 1 - x[::-1]]), lw=4, c='black', alpha=0.7)
    ax.plot(np.concatenate([y, 1 - x[::-1]]), np.concatenate([x, 1 - y[::-1]]), lw=4, c='black', alpha=0.7,
            label='Poisson noise, n = '+str(n))

    ax.legend(loc=2)

    if title:
        ax.set_title(title)

    if savefig:
        fig_filename = savefig
        fig_folder = os.path.dirname(fig_filename)

        mkdirs(fig_folder)
        fig.savefig(fig_filename)
        plt.close(fig)

    else:
        plt.tight_layout()
        plt.ion()
        plt.show()


def plot_minor_allele_example(data, title='', VERBOSE=0, savefig=False):
    '''Plot minor allele in a typical sample'''
    plt.ioff()

    if VERBOSE:
        print 'Plot minor alleles of example sample'

    fig, axs = plt.subplots(1, 2, figsize=(8, 5), sharey=True,
                            gridspec_kw={'width_ratios': [3, 1]})
    sns.set_style('whitegrid')

    labels = ['control', 'patient']
    alphas = [0.6, 1]
    colors = [sns.color_palette()[i] for i in [1, 0]]
    shapes = ['s', 'o']

    for idat, datum in enumerate(data):
        y = datum['freq_minor']
        x = np.arange(len(y))
        #axs[0].plot(x, y, lw=1.5, alpha=0.8)
        axs[0].scatter(x, y,
                       marker=shapes[idat],
                       lw=1.5, edgecolor='none',
                       facecolor=colors[idat],
                       zorder=idat+1)

        h = np.histogram(y, bins=np.logspace(-4, 0, 27))
        axs[1].barh(h[1][:-1], h[0], (h[1][1:] - h[1][:-1]),
                    color=colors[idat],
                    alpha=alphas[idat],
                    zorder=2 - idat)

    axs[0].set_xlabel('Position [bp]')
    axs[0].set_ylabel('Minor allele frequency')
    axs[0].set_yscale('log')
    axs[0].set_ylim(10**(-4), 1)
    axs[0].set_xlim(-20, y.nonzero()[0][-1] + 21)
    axs[0].grid(True)

    axs[1].set_xlabel('Number of alleles')
    axs[1].grid(True)
    axs[1].set_yscale('log')
    axs[1].set_xlim(0.8, 2 * h[0].max())
    axs[1].set_xscale('log')
    #axs[1].set_xticks([0, 50, 100, 150])

    if title:
        fig.suptitle(title)

    if savefig:
        fig_filename = savefig
        fig_folder = os.path.dirname(fig_filename)

        mkdirs(fig_folder)
        fig.savefig(fig_filename)
        plt.close(fig)

    else:
        plt.tight_layout()
        plt.ion()
        plt.show()


def plot_haplotype_tree_example(data, title='', VERBOSE=0, savefig=False,
                                color_by='DSI'):
    '''Plot tree of minor haplotypes in a typical sample'''
    from Bio import Phylo

    def assign_color(tree, cmap='jet', attrname='DSI'):
        '''Assign color to leaves based on a colormap'''
        if isinstance(cmap, basestring):
            from matplotlib import cm
            cmap = getattr(cm, cmap)

        attr_max = max(getattr(leaf, attrname) for leaf in tree.get_terminals())
        def get_color(node):
            return map(int, np.array(cmap(getattr(node, attrname)/attr_max)[:-1]) * 255)
    
        for node in tree.get_terminals():
            node.color = get_color(node)
    
        # For internal nodes, set the attribute (e.g. age) as the arithmetic mean of
        # the children clades
        for node in tree.get_nonterminals(order='postorder'):
            setattr(node, attrname, np.mean([getattr(c, attrname) for c in node.clades]))
            node.color = get_color(node)


    plt.ioff()

    if VERBOSE:
        print 'Plot haplotype tree of example sample'

    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    sns.set_style('white')
    ax.grid(False)

    for datum in data:
        tree = datum['tree']

        times = sorted(set([leaf.DSI for leaf in tree.get_terminals()]))

        assign_color(tree, attrname=color_by)
        labels = ['{:>8s}'.format('{:>3.0%}'.format(leaf.frequency)+' '+
                   '{:>2d}'.format(int(float(leaf.DSI) / 30.5))+
                  'm')
                  for leaf in tree.get_terminals()]
        depths = tree.depths()
        maxdepth = max(depths.itervalues())

        # Collect data for circle plot
        rmin = 5
        rmax = 150
        rfun = lambda hf: rmin + (rmax - rmin) * (hf - 0.015)**(0.5)
        data_circles = []
        for il, leaf in enumerate(tree.get_terminals(), 1):
            it = times.index(leaf.DSI)
            hf = leaf.frequency
            r = rfun(hf)
            y = il
            x = depths[leaf]
            c = [tone / 255.0 for tone in leaf.color.to_rgb()]
            data_circles.append((x, y, 2 * r, c))

        # Draw the tree
        Phylo.draw(tree, show_confidence=False, label_func=lambda x: '', axes=ax,
                   do_show=False)
        ax.set_ylim((ax.get_ylim()[0] + 2, -2))
        ax.set_ylabel('')
        ax.set_yticklabels([])
        for item in ax.get_xticklabels():
            item.set_fontsize(16)
        ax.set_xlabel('Genetic distance [changes / site]', fontsize=16, labelpad=10)

        # Add circles to the leaves
        (x, y, s, c) = zip(*data_circles)
        ax.scatter(x, y, s=s, c=c, edgecolor='none', zorder=2)
        ax.set_xlim(-0.04 * maxdepth, 1.04 * maxdepth)

        # Draw a "legend" for sizes
        datal = [{'hf': 0.05, 'label': '5%'},
                 {'hf': 0.20, 'label': '20%'},
                 {'hf': 1.00, 'label': '100%'}]
        ax.text(0.98 * maxdepth, 1.3, 'Haplotype frequency:', fontsize=16, ha='right')
        for idl, datuml in enumerate(datal):
            r = rfun(datuml['hf'])
            y = 5 + 3 * idl
            ax.scatter(0.85 * maxdepth, y, s=r,
                       facecolor='k',
                       edgecolor='none')
            ax.text(0.98 * maxdepth, y + 0.9, datuml['label'], ha='right')

        # Draw legend for times
        ytext = 0.43 * ax.get_ylim()[0]
        ax.text(0.01 * maxdepth, ytext, 'Time:', fontsize=16)
        datal = [{'time': t, 'color': cm.jet(1.0 * t / max(times))} for t in times]
        for ico, datuml in enumerate(datal):
            y = ytext + 2.8 + 3 * ico
            ax.scatter(0.00 * maxdepth, y, s=rfun(0.5),
                       facecolor=datuml['color'],
                       edgecolor='none')
            ax.text(0.19 * maxdepth, y + 0.9,
                    str(int(datuml['time'] / 30.5))+' months',
                    ha='right')

    plt.tight_layout(rect=(0, 0, 0.98, 1))

    if title:
        fig.suptitle(title)

    if savefig:
        fig_filename = savefig
        fig_folder = os.path.dirname(fig_filename)

        mkdirs(fig_folder)
        fig.savefig(fig_filename)
        plt.close(fig)

    else:
        plt.tight_layout(rect=(0, 0, 0.98, 1))
        plt.ion()
        plt.show()


def plot_substitution_rate(data,
                           regions,
                           title='',
                           VERBOSE=0, savefig=False):
    '''Plot the substitution rates'''

    if VERBOSE:
        print 'Plot substitution rates'
    
    # Sort patient codes
    pcodes = sorted(set(data['pcode']), key=lambda x: int(x[1:]))

    xfun = lambda region, pcode: regions.index(region) + 0.05 * pcodes.index(pcode)
    cfun = lambda pcode: cm.jet(1.0 * pcodes.index(pcode) / len(pcodes))

    fig, ax = plt.subplots(figsize=(0.75 * len(regions), 5))

    for pcode in pcodes:
        datum = (data
                 .loc[data['pcode'] == pcode]
                 .set_index('region', drop=False)
                 .loc[regions])

        y = datum['rate']
        x = map(lambda x: xfun(*x[1]), datum[['region', 'pcode']].iterrows())
        color = cfun(pcode)

        ax.scatter(x, y, color=color, s=90, label=pcode)
        if len(regions) > 1:
            ax.plot(x, y, color=color, lw=1, alpha=0.4)

    ylim = (0.5 * data['rate'].min(), 1.5 * data['rate'].max())
    xticksminor = [xfun(region, pcodes[len(pcodes) //2]) for region in regions]
    xticksmajor = [xt - 0.5 for xt in xticksminor] + [xticksminor[-1] + 0.5]
    xticklabels = regions
    xlim = (xticksminor[0] - 0.04 * (xticksminor[-1] - xticksminor[0]),
            xticksminor[0] + 1.04 * (xticksminor[-1] - xticksminor[0]))


    ax.set_ylim(*ylim)
    ax.set_xlim(*xlim)
    ax.set_xticks(xticksmajor)
    ax.set_xticklabels([])
    ax.set_xticks(xticksminor, minor=True)
    ax.set_xticklabels(xticklabels, fontsize=16, minor=True)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=16)
    ax.grid(True)
    ax.set_yscale('log')

    ax.set_xlabel('Genomic region', fontsize=16)
    ax.set_ylabel('Substitution rate\n[changes / year / site]', labelpad=10, fontsize=16)
    ax.legend(loc='lower left',
              title='Patients',
              ncol=4,
              fontsize=16,
              bbox_to_anchor=(0.10, 0.00))

    if title:
        fig.suptitle(title)

    if savefig:
        fig_filename = savefig
        fig_folder = os.path.dirname(fig_filename)

        mkdirs(fig_folder)
        plt.tight_layout(rect=(0, 0, 0.98, 1))
        fig.savefig(fig_filename)
        plt.close(fig)

    else:
        plt.tight_layout(rect=(0, 0, 0.98, 1))
        plt.ion()
        plt.show()


def plot_substitution_rate_sliding(data,
                                   regions,
                                   title='',
                                   VERBOSE=0, savefig=False):
    '''Plot the substitution rates'''
    from ..analysis.substitution_rate.rate_sliding_window import plot_region_boxes

    if VERBOSE:
        print 'Plot substitution rates in a sliding window'
    
    # Sort patient codes
    pcodes = sorted(set(data['pcode']), key=lambda x: int(x[1:]))

    cfun = lambda pcode: cm.jet(1.0 * pcodes.index(pcode) / len(pcodes))

    fig, axs = plt.subplots(2, 1,
                            sharex=True,
                            figsize=(10, 7),
                            gridspec_kw={'height_ratios':[4, 1]})

    ax = axs[0]
    for pcode in pcodes:
        datum = data.loc[data['pcode'] == pcode].iloc[0]

        y = datum['rate']
        x = datum['x']
        color = cfun(pcode)

        ax.plot(x, y, color=color, lw=2, label=pcode)

    ax.set_ylim(1e-4, 2e-1)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=16)
    ax.grid(True)
    ax.set_yscale('log')
    for item in ax.get_yticklabels():
        item.set_fontsize(16)

    ax.set_ylabel('Substitution rate\n[changes / year / site]', labelpad=20,
                  fontsize=16)
    ax.legend(loc='upper center',
              title='Patients',
              ncol=4,
              fontsize=16)

    ax = axs[1]
    ax.set_ylim(-4, 26)
    ax.set_xlabel('Position in HXB2 [bp]', fontsize=16, labelpad=20)
    ax.set_yticks([])
    for item in ax.get_xticklabels():
        item.set_fontsize(16)
    ax.grid(True, axis='x')

    plot_region_boxes(regions, ax, VERBOSE=VERBOSE)

    if title:
        fig.suptitle(title)

    if savefig:
        fig_filename = savefig
        fig_folder = os.path.dirname(fig_filename)

        mkdirs(fig_folder)
        plt.tight_layout(rect=(0, 0, 0.98, 1), h_pad=0.001)
        fig.savefig(fig_filename)
        plt.close(fig)

    else:
        plt.tight_layout(rect=(0, 0, 0.98, 1), h_pad=0.001)
        plt.ion()
        plt.show()


def plot_mutation_rate(data,
                       title='', savefig=False,
                       VERBOSE=0):
    '''Plot mutation rate estimate'''
    plt.ioff()

    if VERBOSE:
        print 'Plot mutation rate'

    sns.set_style('dark')

    from ..utils.sequence import alphal

    fits = data['fits']
    comp = data['comp']


    # Supplementary figure for the accumulation
    if 'data' in data:
        data = data['data']
        muts = list(np.unique(data.loc[:, 'mut']))
        Sbins = [(0, 0.1), (0.1, 2)]
        Sbin = Sbins[-1]
        datap = (data
                 .loc[(data.loc[:, 'S'] > Sbin[0]) & (data.loc[:, 'S'] <= Sbin[1])]
                 .loc[:, ['mut', 'time', 'af']]
                 .groupby(['mut', 'time'])
                 .mean()
                 .loc[:, 'af']
                 .unstack('time'))

        fig1, ax1 = plt.subplots()
        for mut, datum in datap.iterrows():
            x = np.array(datum.index) / 30.5
            y = np.array(datum)

            # Get rid of masked stuff (this happens if we miss data in the means)
            ind = -(np.isnan(x) | np.isnan(y))
            x = x[ind]
            y = y[ind]

            color = cm.jet(1.0 * muts.index(mut) / len(muts))
            ax1.scatter(x, y, color=color, label=mut)

            m = fits.loc[mut] * 30.5
            xfit = np.linspace(0, 100)
            ax1.plot(xfit, m * xfit, color=color, lw=1.5, alpha=0.5)

        ax1.set_xlabel('Time from infection [months]')
        ax1.set_ylabel('Allele frequency')
        ax1.set_ylim(1e-4, 1)
        ax1.set_yscale('log')
        ax1.grid(True)
        ax1.legend(loc='upper left', fontsize=10, ncol=2)

        plt.tight_layout(rect=(0, 0, 0.98, 1))

    else:
        fig1 = None

    # Plot the mutation rate matrix
    fig2, ax2 = plt.subplots(figsize=(5.3, 4))

    mu = np.ma.masked_all((4, 4))
    for mut, fit in fits.iteritems():
        n1 = mut[0]
        n2 = mut[-1]
        mu[alphal.index(n1), alphal.index(n2)] = fit

    # Plot the log10 for better dynamic range
    vmin = -9
    vmax = -4
    h = ax2.imshow(np.log10(mu),
                   interpolation='nearest',
                   cmap=cm.jet,
                   vmin=vmin, vmax=vmax)
    
    ax2.set_xticks(np.arange(4) + 0.0)
    ax2.set_yticks(np.arange(4) + 0.0)
    ax2.set_xticklabels(alphal[:4], fontsize=16)
    ax2.set_yticklabels(alphal[:4], fontsize=16)
    ax2.set_xlabel('to', fontsize=16)
    ax2.set_ylabel('from', fontsize=16)


    cticks = np.arange(vmin, vmax+1)

    from matplotlib.ticker import Formatter
    class LogTickFormatter(Formatter):
        def __call__(self, x, pos=None):
            '''Transform the logs into their 10**'''
            return '$10^{{{:1.0f}}}$'.format(x)

    cb = plt.colorbar(h, ax=ax2, ticks=cticks, format=LogTickFormatter())
    cb.set_label('Mutation rate\n[changes / position / day]', rotation=270,
                 labelpad=45, fontsize=16)
    cb.ax.tick_params(labelsize=16) 

    plt.tight_layout()

    # Plot the comparison with Abram 2010
    fig3, ax3 = plt.subplots(figsize=(6, 5))
    xmin_exp = -9
    xmax_exp = -4
    xpl = np.logspace(xmin_exp, xmax_exp, 1000)

    x = np.array(comp.loc[:, 'Abram2010'])
    y = np.array(comp.loc[:, 'new'])
    r = np.corrcoef(np.log10(x), np.log10(y))[0, 1]

    ax3.plot(xpl, xpl, color='grey', lw=2)
    ax3.scatter(x, y,
                s=40,
                c='k',
                label='{:2.0%}'.format(r))
    ax3.set_xlabel('Rate from Abram et al. 2010', fontsize=16)
    ax3.set_ylabel('Rate from longitudinal samples', fontsize=16)
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlim(10**(xmin_exp), 10**(xmax_exp))    
    ax3.set_ylim(10**(xmin_exp), 10**(xmax_exp))    
    ax3.legend(loc='upper left', title='Correlation\ncoefficient:', fontsize=16)
    ax3.grid(True)
    for item in ax3.get_xticklabels():
        item.set_fontsize(16)
    for item in ax3.get_yticklabels():
        item.set_fontsize(16)

    plt.tight_layout()

    if title:
        fig1.suptitle(title)
        fig2.suptitle(title)

    if savefig:
        from itertools import izip
        figs = [fig for fig in (fig1, fig2, fig3) if fig is not None]
        for fig, sf in izip(figs, savefig):
            fig_filename = sf
            fig_folder = os.path.dirname(fig_filename)

            mkdirs(fig_folder)
            fig.savefig(fig_filename)
            plt.close(fig)

    else:
        plt.ion()
        plt.show()


def plot_fitness_cost(fits, title='', VERBOSE=0, savefig=False):
    '''Plot the estimated fitness value, for all regions together'''
    fig, ax = plt.subplots(figsize=(8, 6))
    regions = list(set(fits['region']))
    for (region, fitsreg) in fits.groupby('region'):
        ax.plot(fitsreg['S'], fitsreg['s'], lw=2, label=region,
                color=cm.jet(1.0 * regions.index(region) / len(regions)),
                )

    ax.set_xlabel('Entropy in subtype [bits]', fontsize=16)
    ax.set_ylabel('Fitness cost', fontsize=16)
    ax.set_ylim(1e-3, 1)
    ax.set_xlim(1e-3, 2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    for item in ax.get_xticklabels():
        item.set_fontsize(16)
    for item in ax.get_yticklabels():
        item.set_fontsize(16)
    ax.grid(True, which='both')
    ax.legend(loc='upper right', title='Genomic region:', ncol=1,
              fontsize=14)

    if title:
        ax.set_title(title, fontsize=20)

    plt.tight_layout()

    if savefig:
        fig_filename = savefig
        fig_folder = os.path.dirname(fig_filename)

        mkdirs(fig_folder)
        fig.savefig(fig_filename)
        plt.close(fig)

    else:
        plt.ion()
        plt.show()



def plot_allele_freq_example(data, title='', VERBOSE=0, savefig=False):
    '''Plot the estimated fitness value, for all regions together'''
    fig, axs = plt.subplots(1, 3, figsize=(9, 4))
    sns.set_style('darkgrid')

    datum = data[0]
    ind = np.arange(len(datum['times']))
    ind = [0, len(ind) // 2, ind[-1]]

    icons0 = datum['aft'][0].argmax(axis=0)

    for ii, i in enumerate(ind):
        ax = axs[ii]
        time = datum['times'][i]
        af = datum['aft'][i]

        af_min = []
        for pos, afpos in enumerate(af.T):
            afpos = afpos.copy()
            afpos[icons0[pos]] = 0
            afpos.sort()
            af_min.append(afpos[-1])
        af_min = np.array(af_min)

        x = np.arange(af.shape[1])
        color = cm.jet(np.linspace(0, 1, len(x)))

        ax.scatter(x, af_min, s=100, c=color, edgecolor='none')
        
        ax.set_ylim(1e-2, 1.35)
        ax.set_xlim(-5, len(x) + 5)
        ax.set_yscale('log')
        ax.grid(True)
        ax.set_title(str(int(time / 30.5))+' months', fontsize=18)
        for item in ax.get_xticklabels():
            item.set_fontsize(18)

        if ii == 0:
            for item in ax.get_yticklabels():
                item.set_fontsize(18)
        else:
            ax.set_yticklabels([])

    axs[1].set_xlabel('Position [bp]', fontsize=18, labelpad=20)
    axs[0].set_ylabel('Allele frequency', fontsize=18, labelpad=20)


    plt.tight_layout()

    if savefig:
        fig_filename = savefig
        fig_folder = os.path.dirname(fig_filename)

        mkdirs(fig_folder)
        fig.savefig(fig_filename)
        plt.close(fig)

    else:
        plt.ion()
        plt.show()


def plot_divergence_cons_pop(data, VERBOSE=0, savefig=False,
                             include_all=False):
    '''Plot divergence'''
    from itertools import izip

    sns.set_style('darkgrid')

    if not include_all:
        data = data.loc[data['class'] != 'all']

    # NOTE: gp41 has overlaps with tat/rev (hard to tell syn/nonsyn)
    # NOTE: gp120 has the V loops that are hard to call syn/nonsyn
    reg_groups = [['PR', 'IN', 'p15', 'RT'],
                  ['p17', 'p24', 'p6', 'p7'],
                  ['vif', 'vpu', 'vpr', 'nef'],
                 ]

    pcodes = np.unique(data['patient']).tolist()
    pcodes.sort(key=lambda x: int(x[1:]))

    cls = np.unique(data['class'])

    fig, axg = plt.subplots(len(reg_groups), len(cls),
                            figsize=(2 + 4 * len(cls), 4 * len(reg_groups)),
                            sharey=True, sharex=True)

    fig.text(0.41, 0.035, 'Time [days from infection]', fontsize=16)
    fig.text(0.035, 0.6, 'Divergence [changes per site]', rotation=90, ha='center', fontsize=16)
    fig.text(0.32, 0.967, 'nonsyn', ha='center', fontsize=16)
    fig.text(0.75, 0.967, 'syn', ha='center', fontsize=16)


    for irow, (axs, reg_group) in enumerate(izip(axg, reg_groups)):
        axs[-1].set_ylabel(', '.join(reg_group), rotation=270, labelpad=27,
                           fontsize=16)
        axs[-1].yaxis.set_label_position("right")
        for ax in axs:
            ax.grid(True)
        
        datap = (data
                 .loc[data['region'].isin(reg_group)]
                 .groupby(['ctype', 'class', 'patient']))

        for (ctype, cl, pcode), datum in datap:
            if ctype == 'consensus':
                ls = '-'
            else:
                ls = '--'

            if cl == 'all':
                ax = axs[-3]
            elif cl == 'nonsyn':
                ax = axs[-2]
            else:
                ax = axs[-1]

            if (cl == 'nonsyn') and (ctype == 'consensus'):
                label = pcode
            else:
                label = ''

            datump = datum.groupby('time', as_index=False).mean()
            x = datump['time']
            y = datump['div']

            ax.plot(x, y,
                    ls=ls,
                    lw=2,
                    color=cm.jet(1.0 * pcodes.index(pcode) / len(pcodes)),
                    label=label,
                  )

        if irow == 0:
            axs[-2].legend(loc=2, ncol=2, fontsize=14, title='Patients:')


    for ax in axs:
        for item in ax.get_xticklabels():
            item.set_fontsize(16)

    for axs in axg:
        for item in axs[0].get_yticklabels():
            item.set_fontsize(16)

    plt.tight_layout(rect=(0.05, 0.05, 0.98, 0.97))
        
    if savefig:
        fig_filename = savefig
        fig_folder = os.path.dirname(fig_filename)

        mkdirs(fig_folder)
        fig.savefig(fig_filename)
        plt.close(fig)

    else:
        plt.ion()
        plt.show()
