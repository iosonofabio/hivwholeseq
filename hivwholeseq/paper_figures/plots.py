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
def HIVEVO_colormap(kind='website'):
    from scipy.interpolate import interp1d
    maps = {'website': ["#5097BA", "#60AA9E", "#75B681", "#8EBC66", "#AABD52",
                        "#C4B945", "#D9AD3D", "#E59637", "#E67030", "#DF4327"],
            'alternative': ["#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                            "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A"],
           }
    colors = maps[kind]
    rgb_colors = []
    for c in colors:
        rgb_colors.append([int(c[i:i+2],16) for i in [1,3,5]]+[255])
    tmp =interp1d(np.linspace(0,1,len(colors)), np.array(rgb_colors, dtype = float).T/255.0)
    cmap = lambda x: [c for c in tmp(x)]
    return cmap


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
    if VERBOSE >= 2:
        print 'Plot allele frequency in overlaps'

    sns.set_style('darkgrid')
    colors = sns.color_palette('Set1', 5)
    fs = 16
    xmin = 1e-3
    
    fig, ax = plt.subplots(figsize=(6, 5))

    legend = set()
    for ida, datum in enumerate(data):
        afjoint = datum['af']
        color = colors[datum['io']]
        if datum['overlap'] not in legend:
            label = datum['overlap']
            #legend.add(datum['overlap'])
        else:
            label = None
        ax.scatter(afjoint[0].ravel(), afjoint[1].ravel(),
                   s=50,
                   color=color,
                   alpha=0.7,
                   #label=label,
                   edgecolor='none')

        # Plot stddev in Poisson sampling
        n = datum['n']
        x = np.linspace(np.log10(xmin), 0, 1000)
        x = 1.0 / (1 + 10**(-x))
        y = x - np.sqrt(x / n)
        ax.plot(np.concatenate([x, 1 - y[::-1]]), np.concatenate([y, 1 - x[::-1]]), lw=3, c=color, alpha=0.5)
        ax.plot(np.concatenate([y, 1 - x[::-1]]), np.concatenate([x, 1 - y[::-1]]), lw=3, c=color, alpha=0.5,
                label=datum['overlap'])

    
    ax.plot([xmin, 1 - xmin], [xmin, 1 - xmin], lw=2, color='k', alpha=0.5)

    ax.set_xlabel('SNV frequency leading fragment', fontsize=fs)
    ax.set_ylabel('SNV frequency trailing fragment', fontsize=fs)
    ax.grid(True)

    if use_logit:
        ax.set_xscale('logit')
        ax.set_yscale('logit')
        ax.set_xlim(xmin, 1 - xmin)
        ax.set_ylim(xmin, 1 - xmin)

    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)
    ax.legend(loc=2, fontsize=fs)

    if title:
        ax.set_title(title)

    plt.tight_layout()

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
    colors = [sns.color_palette()[i] for i in [2, 0]]
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

    axs[0].set_xlabel('Position [bp]', fontsize=18)
    axs[0].set_ylabel('SNV frequency', fontsize=18)
    axs[0].set_yscale('log')
    axs[0].set_ylim(10**(-4), 1)
    axs[0].set_xlim(-20, y.nonzero()[0][-1] + 21)
    axs[0].grid(True)
    axs[0].tick_params(axis='both', labelsize=18)

    axs[1].set_xlabel('Number of positions', fontsize=18)
    axs[1].grid(True)
    axs[1].set_yscale('log')
    axs[1].set_xlim(0.8, 2 * h[0].max())
    axs[1].set_xscale('log')
    axs[1].tick_params(axis='x', labelsize=18)
    #axs[1].set_xticks([0, 50, 100, 150])

    plt.tight_layout(rect=(0, 0, 1, 1))

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
                                color_by='DSI', legend=None):
    '''Plot tree of minor haplotypes in a typical sample'''
    from itertools import izip
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

    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    sns.set_style('white')
    colormap = HIVEVO_colormap()
    for datum, ax in izip(data, axs.ravel()):
        ax.grid(False)
        tree = datum['tree']
        region = datum['region']
        tree.root.branch_length = 0.001
        times = sorted(set([leaf.DSI for leaf in tree.get_terminals()]))

        assign_color(tree, attrname=color_by, cmap=colormap)
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
            cs = [tone / 255.0 * 0.7 for tone in leaf.color.to_rgb()]
            data_circles.append((x, y, 2 * r, c, cs))

        # Draw the tree
        Phylo.draw(tree, show_confidence=False, label_func=lambda x: '', axes=ax,
                   do_show=False)
        ax.set_ylim((1.02 * ax.get_ylim()[0], -0.2 - 1.8 * int(ax.get_ylim()[0] > 30)))
        ax.set_ylabel('')
        ax.set_yticklabels([])
        ax.set_axis_off()
        for item in ax.get_xticklabels():
            item.set_fontsize(16)

        if region != 'V3':
            ax.set_xlabel('Genetic distance [changes / site]', fontsize=16, labelpad=10)
        else:
            ax.set_xlabel('')

        # Add circles to the leaves
        (x, y, s, c, cs) = zip(*data_circles)
        ax.scatter(x, y, s=s, c=c, edgecolor=cs, zorder=2, clip_on=False)
        ax.set_xlim(-0.04 * maxdepth, 1.04 * maxdepth)

        if region == legend:
            # Draw a "legend" for sizes
            datal = [{'hf': 0.05, 'label': '5%'},
                     {'hf': 0.20, 'label': '20%'},
                     {'hf': 1.00, 'label': '100%'}]
            ax.text(0.98 * maxdepth, 1.3, 'Haplotype frequency:', fontsize=16, ha='right')
            for idl, datuml in enumerate(datal):
                r = rfun(datuml['hf'])
                y = 5 + 3 * idl
                ax.scatter(0.80 * maxdepth, y, s=r,
                           facecolor='k',
                           edgecolor='none')
                ax.text(0.98 * maxdepth, y + 0.9, datuml['label'], ha='right')

            # Draw legend for times
            datal = [{'time': t,'color': [tone for tone in colormap(1.0 * t / max(times))[:3]],
            'colorstroke': [tone * 0.7 for tone in colormap(1.0 * t / max(times))[:3]]} for t in times]
            xtext = -0.015
            ytext = (0.93 - 0.06 * min(8, len(datal))) * ax.get_ylim()[0]
            ax.text(xtext-0.02*maxdepth, ytext, 'Time:', fontsize=16)
            for ico, datuml in enumerate(datal):
                ytext += 0.06 * (ax.get_ylim()[0] - ax.get_ylim()[1])

                if ico == 8:
                    ytext -= 4 * 0.06 * (ax.get_ylim()[0] - ax.get_ylim()[1])
                    xtext += 0.28 * maxdepth

                ax.scatter(xtext, ytext, s=rfun(0.5),
                           facecolor=datuml['color'],
                           edgecolor=datuml['colorstroke'],
                           clip_on=False)
                ax.text(xtext + 0.26 * maxdepth, ytext + 0.02 * ax.get_ylim()[0],
                        str(int(datuml['time'] / 30.5))+' months',
                        ha='right',
                        fontsize=14,
                       )

        # Draw scale bar
        xbar = (0.35 + 0.3 * (len(times) >= 9))*(legend == region) * maxdepth
        ybar = 0.90 * ax.get_ylim()[0]
        lbar = 0.15 * maxdepth
        lbar_label = '{:.1G}'.format(lbar)
        lbar = float(lbar_label)
        ax.plot([xbar, xbar + lbar], [ybar, ybar], lw=4, c='k')
        ax.text(xbar + 0.5 * lbar, ybar + 0.11 * ax.get_ylim()[0],
                lbar_label, fontsize=18,
                ha='center')
    plt.tight_layout(rect=(0, 0, 0.98, 1), pad=0, w_pad=0, h_pad=0)

    if title:
        fig.suptitle(title)

    if savefig:
        fig_filename = savefig
        fig_folder = os.path.dirname(fig_filename)

        mkdirs(fig_folder)
        fig.savefig(fig_filename)
        plt.close(fig)

    else:
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
    ax.set_ylabel('Divergence rate\nin observed period\n[changes / year / site]', labelpad=10, fontsize=16)
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
    if VERBOSE:
        print 'Plot substitution rates in a sliding window'

    data_sweeps = data['substitutions']
    data_ctl = data['ctl']
    data = data['substitution_rate']

    # Plot parameters
    fs = 16
    
    def plot_region_boxes(regions, ax, refname='HXB2', VERBOSE=0):
        '''Plot boxes for genomic regions of HXB2'''
        import pandas as pd
        from matplotlib.patches import Rectangle
    
        from hivwholeseq.reference import load_custom_reference
        refseq = load_custom_reference(refname, format='gb')
    
        data = []
        
        y1 = 5
        height = 5
        pad = 2
        for feature in refseq.features:
            if feature.id in regions:
                x = [max(547, feature.location.nofuzzy_start),
                     min(9591, feature.location.nofuzzy_end)]
                name = feature.id
                if name == 'RRE':
                    name = '         RRE'
                elif name == 'V5':
                    name = 'V5   '

                data.append({'name': name,
                             'x1': x[0], 'x2': x[1], 'width': x[1] - x[0]})
    
        data = pd.DataFrame(data)
        data.sort('x1', inplace=True)
        data.index = np.arange(len(data))
        data['height'] = height
        data['parity'] = ((1 + np.arange(len(data))) % 2)
        data['row'] = 'bottom'; data.loc[np.array(data['parity'] == 1), 'row'] = 'top'
        data['y1'] = y1 + (height + pad) * data['parity']
        data['y2'] = data['y1'] + data['height']
    
        for _, datum in data.iterrows():
            r = Rectangle((datum['x1'], datum['y1']),
                          datum['width'],
                          datum['height'],
                          facecolor=[0.8] * 3,
                          edgecolor='k',
                          label=datum['name'])
            
            xt = datum['x1'] + 0.5 * datum['width']
            yt = datum['y1'] + 0.5 * datum['height']
            if datum['row'] == 'top':
                yt +=  height + 0
            else:
                yt -= height + 4.2
    
            ax.add_patch(r)
            ax.text(xt, yt,
                    datum['name'],
                    color='k', 
                    fontsize=fs,
                    ha='center')

            if 'RRE' in datum['name']:
                ax.plot([datum['x1'] + 0.5 * datum['width'],
                         datum['x1'] + 2.5 * datum['width'],
                         datum['x1'] + 2.5 * datum['width'],
                        ],
                        [datum['y1'] + 0.5 * datum['height'],
                         datum['y1'] + 0.5 * datum['height'],
                         datum['y1'] + 0.0 * datum['height']
                        ],
                        lw=2, color='k')


    # Sort patient codes
    pcodes = sorted(set(data['pcode']), key=lambda x: int(x[1:]))
    cfun = lambda pcode: cm.jet(1.0 * pcodes.index(pcode) / len(pcodes))

    fig, axs = plt.subplots(3, 1,
                            sharex=True,
                            figsize=(8, 8),
                            gridspec_kw={'height_ratios':[4, 1.5, 1]})

    # Substitution rate
    ax = axs[0]
    for pcode in pcodes:
        datum = data.loc[data['pcode'] == pcode].iloc[0]

        y = datum['rate']
        x = datum['x']
        color = cfun(pcode)

        ax.plot(x, y, color=color, lw=2, label=pcode)

    ax.set_ylim(1e-4, 2e-1)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=fs)
    ax.grid(True)
    ax.set_yscale('log')
    ax.xaxis.set_tick_params(labelsize=fs)

    ax.set_ylabel('Evolutionary rate\nduring observed period\n[changes / year / site]',
                  labelpad=15,
                  fontsize=fs)
    ax.legend(loc='upper center',
              bbox_to_anchor=(0.5, 0.93),
              ncol=4,
              fontsize=fs)
    ax.text(0.5, 0.93, 'Patients:', fontsize=18,
            ha='center',
            transform=ax.transAxes,
           )

    
    # Substitutions
    ax = axs[1]
    pnames = data_sweeps['pcode'].unique().tolist()
    Lp = len(pnames)

    for pname, datum in data_sweeps.groupby('pcode'):
        x = np.array(datum['pos_ref'])
        y = np.repeat(pnames.index(pname), len(x))

        for ind, marker, s in [(datum['epitope'], 'o', 30),
                               (-datum['epitope'], 'x', 30)]:
            ind = np.array(ind)
            ax.scatter(x[ind], y[ind],
                       s=s, marker=marker,
                       color=cfun(pname),
                      )

    for pname, datum in data_ctl.groupby('pcode'):
        y = pnames.index(pname) + 0.35
        for _, datump in datum.iterrows():
            x_left = datump['start_HXB2']
            x_right = datump['end_HXB2']
            width = x_right - x_left
            ax.plot([x_left, x_right], [y] * 2,
                    color=cfun(pname),
                    lw=2,
                   )

    ax.set_xlim(-50, data_sweeps['pos_ref'].max() + 200)
    ax.set_ylim(Lp, -1)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel('Substitutions', fontsize=fs, labelpad=60)
    ax.grid(True)


    # Genome map
    ax = axs[2]
    ax.set_ylim(-4, 26)
    ax.set_xlabel('Position [bp]', fontsize=fs, labelpad=20)
    ax.set_yticks([])
    ax.xaxis.set_tick_params(labelsize=fs)
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
    fig, axs = plt.subplots(2, 3, figsize=(7, 5))
    sns.set_style('darkgrid')
    fs = 18
    datum = data[0]
    ind = np.arange(len(datum['times']))
    ind = [0, len(ind) // 2, ind[-1]]

    icons0 = datum['aft'][0].argmax(axis=0)

    cmap = HIVEVO_colormap(kind='alternative')
    x = np.arange(datum['aft'].shape[2])
    color = [[float(tmp) for tmp in cmap(p)] for p in np.linspace(0, 1, len(x))]
    for ii, i in enumerate(ind):
        ax = axs[0][ii]
        time = datum['times'][i]
        af = datum['aft'][i]

        af_min = []
        for pos, afpos in enumerate(af.T):
            afpos = afpos.copy()
            afpos[icons0[pos]] = 0
            afpos.sort()
            af_min.append(afpos[-1])
        af_min = np.array(af_min)
        ax.scatter(x, af_min, s=100, c=color, edgecolor='none')
        
        ax.set_ylim(1e-2, 1.35)
        ax.set_xlim(-5, len(x) + 5)
        ax.set_xticks(range(0, len(x), 150))
        ax.set_yscale('log')
        ax.grid(True)
        ax.set_title(str(int(time / 30.5))+' months', fontsize=fs)
        for item in ax.get_xticklabels():
            item.set_fontsize(fs)

        if ii == 0:
            for item in ax.get_yticklabels():
                item.set_fontsize(fs)
        else:
            ax.set_yticklabels([])

    axs[0][1].set_xlabel('Position [bp]', fontsize=fs, labelpad=5)
    fig.text(0.035, 0.5, 'SNV frequency', ha='center', va='center', rotation='vertical',
             fontsize=fs)

    ax = plt.subplot2grid((2, 3), (1, 0), colspan=3)
    tmonth = datum['times']/30.5
    for pos in xrange(datum['aft'].shape[2]):
        for nuc in xrange(4):
            traj = datum['aft'][:,nuc,pos]
            traj[traj<0.003] = 0.003
            if (traj[0] < 0.5) and (traj.max() > 0.05):
                ax.plot(tmonth, traj, c=color[pos])

    ax.set_ylim(1e-2, 1.35)
    ax.set_xlim(0, tmonth[-1] + 1)
    ax.set_yscale('log')
    ax.set_xlabel('Time since infection [months]', fontsize=fs)
    for item in ax.get_xticklabels() + ax.get_yticklabels():
        item.set_fontsize(fs)

    plt.tight_layout(rect=(0.07, 0.02, 0.98, 0.98), pad=0.05, h_pad=0.5, w_pad=0.4)
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

    fig.text(0.41, 0.035, 'Time [months from infection]', fontsize=16)
    fig.text(0.035, 0.6, 'Divergence [changes per site]', rotation=90, ha='center', fontsize=16)

    if len(cls) == 2:
        fig.text(0.32, 0.967, 'nonsyn', ha='center', fontsize=16)
        fig.text(0.75, 0.967, 'syn', ha='center', fontsize=16)
    else:
        fig.text(0.22, 0.967, 'all', ha='center', fontsize=16)
        fig.text(0.53, 0.967, 'nonsyn', ha='center', fontsize=16)
        fig.text(0.82, 0.967, 'syn', ha='center', fontsize=16)
 

    # Bin data in time
    from ..utils.pandas import add_binned_column
    _, times = add_binned_column(data, 'tbin', 'time',
                                 bins=np.linspace(0, 3300, 8),
                                 clip=True)
    data.loc[:, 'tbinned'] = times[data['tbin']]

    for irow, (axs, reg_group) in enumerate(izip(axg, reg_groups)):
        axs[-1].set_ylabel(', '.join(reg_group), rotation=270, labelpad=27,
                           fontsize=16)
        axs[-1].yaxis.set_label_position("right")
        for ax in axs:
            ax.grid(True)
        
        def plot_single(ctype, cl, pcode, datum, groupby='time'):
            '''Plot single curve'''
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

            if pcode in pcodes:
                color = cm.jet(1.0 * pcodes.index(pcode) / len(pcodes))
                alpha = 0.4
            else:
                color = 'k'
                alpha = 1.0

            datump = datum.groupby(groupby, as_index=False).mean()
            x = datump[groupby] / 30.5
            y = datump['div']

            ax.plot(x, y,
                    ls=ls,
                    lw=2,
                    alpha=alpha,
                    color=color,
                    label=label,
                  )


        # Plot single patients
        datap = (data
                 .loc[data['region'].isin(reg_group)]
                 .groupby(['ctype', 'class', 'patient']))
        for (ctype, cl, pcode), datum in datap:
            plot_single(ctype, cl, pcode, datum)

        # Plot average over patients
        datap = (data
                 .loc[data['region'].isin(reg_group)]
                 .groupby(['ctype', 'class']))
        for (ctype, cl), datum in datap:
            plot_single(ctype, cl, 'avg', datum, groupby='tbinned')

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


def plot_divergence_cons_pop_diversity(data, VERBOSE=0, savefig=False,
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

    fig.text(0.41, 0.035, 'Time [months from infection]', fontsize=16)
    fig.text(0.035, 0.6, 'Divergence [changes per site]', rotation=90, ha='center', fontsize=16)

    if len(cls) == 2:
        fig.text(0.32, 0.967, 'nonsyn', ha='center', fontsize=16)
        fig.text(0.75, 0.967, 'syn', ha='center', fontsize=16)
    else:
        fig.text(0.22, 0.967, 'all', ha='center', fontsize=16)
        fig.text(0.53, 0.967, 'nonsyn', ha='center', fontsize=16)
        fig.text(0.82, 0.967, 'syn', ha='center', fontsize=16)

    # Bin data in time
    from ..utils.pandas import add_binned_column
    _, times = add_binned_column(data, 'tbin', 'time',
                                 bins=np.linspace(0, 3300, 8),
                                 clip=True)
    data.loc[:, 'tbinned'] = times[data['tbin']]

    for irow, (axs, reg_group) in enumerate(izip(axg, reg_groups)):
        axs[-1].set_ylabel(', '.join(reg_group), rotation=270, labelpad=27,
                           fontsize=16)
        axs[-1].yaxis.set_label_position("right")
        for ax in axs:
            ax.grid(True)
        
        def plot_single(obs, ctype, cl, pcode, datum, groupby='time'):
            '''Plot single curve'''
            if obs == 'diversity':
                ls = '-'
                dashes = [8, 4, 2, 4, 2, 4]
            elif ctype == 'consensus':
                ls = '-'
                dashes = []
            else:
                ls = '--'
                dashes = [8, 6]

            if cl == 'all':
                ax = axs[-3]
            elif cl == 'nonsyn':
                ax = axs[-2]
            else:
                ax = axs[-1]

            if (cl == 'nonsyn') and (ctype == 'consensus'):
                label = pcode
            elif (cl == 'syn') and (pcode == 'avg'):
                if obs == 'diversity':
                    label = 'Diversity'
                elif ctype == 'consensus':
                    label = 'Divergence (consensus)'
                else:
                    label = 'Divergence (population)'
            else:
                label = ''

            if pcode in pcodes:
                color = cm.jet(1.0 * pcodes.index(pcode) / len(pcodes))
                alpha = 0.25
                lw = 1
            else:
                color = 'k'
                alpha = 1.0
                lw = 2

            datump = datum.groupby(groupby, as_index=False).mean()
            x = datump[groupby] / 30.5
            y = datump['div']

            ax.plot(x, y,
                    ls=ls,
                    dashes=dashes,
                    lw=lw,
                    alpha=alpha,
                    color=color,
                    label=label,
                  )


        # Plot single patients
        datap = (data
                 .loc[data['region'].isin(reg_group)]
                 .groupby(['obs', 'ctype', 'class', 'patient']))
        for (obs, ctype, cl, pcode), datum in datap:
            plot_single(obs, ctype, cl, pcode, datum)

        # Plot average over patients
        datap = (data
                 .loc[data['region'].isin(reg_group)]
                 .groupby(['obs', 'ctype', 'class']))
        for (obs, ctype, cl), datum in datap:
            plot_single(obs, ctype, cl, 'avg', datum, groupby='tbinned')

        if irow == 0:
            axs[-2].legend(loc=2, ncol=2, fontsize=14, title='Patients:')
            axs[-1].legend(loc=2, ncol=1, fontsize=14)


    for ax in axs:
        ax.xaxis.set_tick_params(labelsize=16)

    for axs in axg:
        ax.yaxis.set_tick_params(labelsize=16)

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


def plot_sfs_syn_nonsyn(data, VERBOSE=0, savefig=False):
    '''Plot SFS synonymous/nonsynonymous'''
    if VERBOSE >= 2:
        print 'Plot SFS syn/nonsyn'

    sns.set_style('darkgrid')
    fs = 18

    plot_props = [{'mutclass': 'syn',
                   'color': 'steelblue',
                   'label': 'synonymous',
                  },
                  {'mutclass': 'nonsyn',
                   'color': 'darkred',                      
                   'label': 'nonsynonymous',
                  }]

    fig, ax = plt.subplots(figsize=(5, 4)) 
    for ip, props in enumerate(plot_props):
        mutclass = props['mutclass']
        color = props['color']
        label = props['label']

        arr = data.loc[:, mutclass+' fraction']
        y = np.array(arr)
        x = np.array(data.loc[:, 'afbin_left'])
        w = np.array(data.loc[:, 'afbin_right'] - x)

        n_props = len(plot_props)
        ws = 1.0 / (n_props + 1)
        x += ws * w * ip

        bottom = 1e-3
        height = y - bottom
        ax.bar(x, height, width=ws * w, bottom=bottom,
               color=color,
               label=label)


    ax.set_xlabel('Allele frequency', fontsize=fs)
    ax.set_xlim(0, 1)
    ax.set_yscale('log')
    ax.grid(True)
    ax.legend(loc='upper right', ncol=1, fontsize=fs)
    ax.set_ylabel('Fraction of of SNVs', fontsize=fs)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)

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


def plot_divergence_diversity(data, VERBOSE=0, savefig=False):
    '''Plot divergence and diversity in various genomic regions'''
    from itertools import izip

    if VERBOSE >= 2:
        print 'Plot divergence and diversity'

    colors = ["#5097BA", "#60AA9E", "#75B681", "#8EBC66", "#AABD52", "#C4B945",
              "#D9AD3D", "#E59637", "#E67030", "#DF4327"]

    # NOTE: gp41 has overlaps with tat/rev (hard to tell syn/nonsyn)
    # NOTE: gp120 has the V loops that are hard to call syn/nonsyn
    reg_groups = [['PR', 'IN', 'p15', 'RT'],
                  ['p17', 'p24', 'p6', 'p7'],
                  ['vif', 'vpu', 'vpr', 'nef'],
                  ['gp120_noVloops'],
                 ]
    fs = 17
    alpha = 0.8
    pcodes = np.unique(data['patient']).tolist()
    pcodes.sort(key=lambda x: int(x[1:]))
    cls = np.unique(data['class'])
    fig, axs = plt.subplots(1, 2,
                            figsize=(8, 4),
                            sharey=True,
                            sharex=True)

    fig.text(0.38, 0.035, 'Time [days from infection]', fontsize=fs)

    axs[0].set_title('nonsyn', ha='center', fontsize=fs)
    axs[1].set_title('syn', ha='center', fontsize=fs)

    # Bin data in time
    from hivwholeseq.utils.pandas import add_binned_column
    _, times = add_binned_column(data, 'tbin', 'time',
                                 bins=np.linspace(0, 3300, 8),
                                 clip=True)
    data.loc[:, 'tbinned'] = times[data['tbin']]
    region_data = [(data
                    .loc[data['region'].isin(reg)]
                    .groupby(['obs', 'ctype', 'class']))
                    for reg in reg_groups]


    def label_func(dclass, obs, reg_label):
        if dclass=='syn' and reg_label=='accessory':
            return obs
        elif dclass=='nonsyn' and obs=='divergence':
            return reg_label
        else:
            return None


    for ax, dclass, ls2 in izip(axs, ['nonsyn', 'syn'], ['-', '--']):
        ax.grid(True)
        for ctype, obs, ls in izip(['population', 'population'], ['divergence','diversity'], ['-','--']):
            for reg, reg_label, color in izip(region_data, 
                                    ['enzymes','structural', 'accessory', 'envelope'], 
                                    [colors[i] for i in [0,4,9,7]]):
                tmp = reg.get_group((obs, ctype, dclass)).groupby('tbinned', as_index=False).mean()
                x = tmp['tbinned']
                y = tmp['div']
                ax.plot(x, y,
                    ls=ls,
                    lw=2,
                    color=color,
                    alpha=alpha,
                    label=label_func(dclass, obs, reg_label),
                  )
        
        ax.legend(loc=2, fontsize=fs)

    ax.set_ylim([0,0.04])
    axs[0].set_yticks([0,0.025, 0.05])
    axs[1].set_yticks([0,0.025, 0.05])
    axs[0].set_xticks([0, 1500, 3000])
    axs[1].set_xticks([0, 1500, 3000])
    axs[0].yaxis.set_tick_params(labelsize=fs)
    axs[0].xaxis.set_tick_params(labelsize=fs)
    axs[1].xaxis.set_tick_params(labelsize=fs)
    #axs[0].set_ylabel('Divergence', fontsize=fs, labelpad=10)
    #axs[1].yaxis.set_label_position('right')
    #axs[1].set_ylabel('Diversity', fontsize=fs, labelpad=10)
    plt.tight_layout(rect=(0, 0.05, 0.98, 1))

    if savefig:
        fig_filename = savefig
        fig_folder = os.path.dirname(fig_filename)

        mkdirs(fig_folder)
        fig.savefig(fig_filename)
        plt.close(fig)

    else:
        plt.ion()
        plt.show()



def plot_LD(data, VERBOSE=0, savefig=False):
    '''Plot linkage disequilibrium'''
    if VERBOSE >= 2:
        print 'Plot LD'

    sns.set_style('darkgrid')
    fs=16

    binc = data['binc']
    Dp_vs_distance = data['Dp']

    fig = plt.figure(2, figsize=(5.5, 4.3))
    ax = plt.subplot(111)
    for frag, y in Dp_vs_distance.iteritems():
        ax.plot(binc, y, label=frag, lw=3)

    ax.set_xticks(range(0, 401, 100))
    ax.set_yticks(np.arange(0., 1.01, 0.2))
    for item in ax.get_xticklabels() + ax.get_yticklabels():
        item.set_fontsize(fs)
    ax.set_ylabel("linkage disequilibrium D'", fontsize=fs)
    ax.set_xlabel("distance [bp]", fontsize=fs)
    ax.legend(loc=1, fontsize=fs, ncol=3)

    plt.tight_layout(rect=(0, 0, 0.98, 1))

    if savefig:
        fig_filename = savefig
        fig_folder = os.path.dirname(fig_filename)

        mkdirs(fig_folder)
        fig.savefig(fig_filename)
        plt.close(fig)

    else:
        plt.ion()
        plt.show()


def plot_coverage_example(data, VERBOSE=0, savefig=False):
    '''Plot an example of coverage'''
    if VERBOSE >= 2:
        print 'Plot coverage'

    sns.set_style('darkgrid')
    fs=16
    colormap = HIVEVO_colormap(kind='alternative')

    fig, ax = plt.subplots(figsize=(8, 4))
    L = len(data['cov']['genomewide']['cov'])

    y = data['n_templates']
    ax.plot([0, L], [y] * 2, lw=10, color='darkred',
            label='Number of HIV-1 RNA molecules in serum sample',
            zorder=1,
           )

    x = np.arange(L)
    ax.plot(x, data['cov']['genomewide']['cov'],
            color='k',
            lw=2,
            zorder=10,
            )

    for ifr, fragment in enumerate(['F'+str(i) for i in xrange(1, 7)]):
        y = data['cov'][fragment]['cov']
        x = np.arange(len(y)) + data['cov'][fragment]['pos']
        ax.scatter(x, y, color=colormap(1.0 * ifr / 6),
                   s=65,
                   marker='o',
                   zorder=3,
                  )

    ax.legend(loc='lower right', fontsize=fs)
    ax.set_xlabel('Position [bp]', fontsize=fs)
    ax.set_ylabel('Coverage', fontsize=fs)
    ax.set_yscale('log')
    ax.set_ylim(0.105, 2e5)
    ax.set_xlim(-200, L + 200)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)

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


def plot_af_above_threshold(datap, threshold=0.01, singles=False,
                            VERBOSE=0,
                            savefig=False):
    '''Plot the fraction of alleles above a certain threshold'''
    if VERBOSE >= 2:
        print 'Plot allele frequencies above thresholds'

    Ssubs = datap.Ssub.unique().tolist()

    sns.set_style('darkgrid')
    colormap = cm.jet
    fs = 14

    fig, ax = plt.subplots(figsize=(5.5, 4))

    if singles:
        for (Ssub, pcode), datump in datap.groupby(['Ssub', 'pcode']):
            x = np.array(datump['time'])
            y = np.array(datump['frac > '+'{:.1G}'.format(threshold)])

            ax.plot(x, y,
                    lw=1,
                    ls='-', marker='x', ms=10,
                    color=colormap(1.0 * Ssubs.index(Ssub) / len(Ssubs)),
                    alpha=0.7,
                   )

    # Average over patients
    for Ssub, datump in datap.groupby('Ssub'):
        datumpt = datump.groupby('time').mean()
        ddatumpt = datump.groupby('time').std()
        x = np.array(datumpt.index)
        y = np.array(datumpt['frac > '+'{:.1G}'.format(threshold)])
        dy = np.array(ddatumpt['frac > '+'{:.1G}'.format(threshold)])
        label = '{:.1G}'.format(Ssub)
        ax.errorbar(x, y, yerr=dy,
                    lw=2,
                    ls='-', marker='o', ms=12,
                    color=colormap(1.0 * Ssubs.index(Ssub) / len(Ssubs)),
                    label=label,
                   )

    # Change legend order
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1],
              loc=2, title='Subtype B entropy:', fontsize=fs)

    ax.set_ylim(ymin=-0.02)
    ax.set_xlabel('Time [days from infection]', fontsize=fs)
    ax.set_ylabel('Fraction of SNVs > '+'{:.1G}'.format(threshold), fontsize=fs)
    ax.grid(True)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)

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


def plot_af_entropy_awayto(data, VERBOSE=0, savefig=False):
    '''Plot average frequency (simpler plot)'''
    bins = data['bins']
    datap = data['af_avg']

    bins_S, binsc_S = bins['entropy']

    fs = 16
    sns.set_style('darkgrid')

    fig, ax = plt.subplots(figsize=(5, 4))
    keys = ('to', 'away')
    for awayto in keys:
        datum = datap.groupby('awayto').get_group(awayto)
        if awayto == 'away':
            ls = '--'
            color = 'goldenrod'
            marker = 's'
        else:
            ls = '-'
            color = 'darkolivegreen'
            marker = 'v'

        x = binsc_S[datum['Sbin']]
        y = datum['af']
        dy = datum['af_std']

        ax.errorbar(x, y, yerr=dy,
                    ls=ls, lw=2,
                marker=marker,
                markersize=10,
                label=awayto, color=color)

    ax.set_xlabel('Entropy in subtype [bits]', fontsize=fs)
    ax.set_ylabel('Average SNV frequency', fontsize=fs)
    ax.set_ylim(1e-4, 1e-1)
    ax.set_xlim(1e-3, 5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='upper left', fontsize=fs)
    ax.grid(True)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)
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


def plot_n_muts_awayto(data, VERBOSE=0, savefig=False):
    '''Plot the averages'''
    if VERBOSE >= 2:
        print 'Plot number of mutations away/to subtype B consensus'

    import pandas as pd

    n_muts = data['n_muts']
    frac = data['fraction']
    frac0 = data['fraction0']

    sns.set_style('darkgrid')
    fs = 16
    colormap = cm.jet
    fig, ax = plt.subplots(figsize=(5, 4))

    # Plot sum
    x = np.array(n_muts['time'])
    y = np.array(n_muts['nmuts'])
    dy = np.array(n_muts['nmuts_std'])
    ax.errorbar(x, y, yerr=dy,
            lw=2,
            marker='o',
            markersize=15,
            color='steelblue',
           )

    ax.set_xlabel('Time from infection [days]', fontsize=fs)
    ax.set_ylabel('N. mutations / genome [no env]', fontsize=fs, labelpad=15)
    ax.set_xticks([0, 1000, 2000, 3000])
    ax.set_yticks([0, 50, 100, 150, 200])
    ax.set_ylim(0, 200)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)
    ax.grid(True)

    # Plot fraction (from ratio)
    ax = ax.twinx()
    x = np.array(frac['time'])
    y = 100. * np.array(frac['fraction to/away'])
    dy = 100. * np.array(frac['fraction_std'])
    ax.errorbar(x, y, yerr=dy,
            lw=2,
            marker='o',
            markersize=15,
            color='darkred',
           )
    y0 = 100. * np.repeat(frac0, len(x))
    ax.plot(x, y0, lw=2, c='k', ls='--')

    yticks = [0, 20, 40]
    ax.set_ylim(yticks[0], yticks[-1])
    ax.set_yticks(yticks)
    ax.set_yticklabels([str(x)+'%' for x in yticks])
    ax.set_xlim(-100, 3100)
    ax.set_ylabel('SNVs to subtype B consensus',
                  fontsize=fs,
                  rotation=90,
                  labelpad=15)
    ax.yaxis.set_tick_params(labelsize=fs)
    ax.grid(True)

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


def plot_entropy_correlation(datap, VERBOSE=0, savefig=False):
    '''Plot the correlation between intrapatient and subtype diversity'''
    if VERBOSE >= 2:
        print 'Correlation in Shannon entropy between patient and subtype'''
    
    sns.set_style('darkgrid')
    colormap = cm.jet#HIVEVO_colormap('website')
    fs = 16

    datump = datap.groupby(['pcode', 'time'], as_index=False).mean()
    datump['rho_std'] = np.array(datap.groupby(['pcode', 'time']).std()['rho'])
    pcodes = datap['pcode'].unique().tolist()
    pcodes.sort(key=lambda x: int(x[1:]))

    fig, ax = plt.subplots(figsize=(6, 5))
    for ip, pcode in enumerate(pcodes):
        datumpp = datump.groupby('pcode').get_group(pcode)
        x = np.array(datumpp['time'])
        y = np.array(datumpp['rho'])
        dy = np.array(datumpp['rho_std'])
        ax.errorbar(x, y, dy,
                    ls="none",
                    marker='o',
                    markersize=15,
                    color=colormap(1.0 * ip / len(pcodes)),
                    label=pcode,
                   )

    # Plot linear fit excluding p7 which is very late and slow
    ind = (datump['pcode'] != 'p7') & (-np.isnan(datump['rho_std']))
    x = datump['time'].loc[ind]
    y = datump['rho'].loc[ind]
    dy = datump['rho_std'].loc[ind]
    m, q = np.polyfit(x, y, 1, w=1.0/dy)
    xfit = np.linspace(0, 3500, 1000)
    yfit = q + m * xfit
    ax.plot(xfit, yfit, lw=2, c='k')

    ax.legend(loc='upper left', fontsize=fs-3, ncol=2, title='Patient:')
    ax.set_xlim(-100, 3500)
    ax.set_ylim(-0.03, 0.70)
    ax.set_xlabel('Time from infection [days]', fontsize=fs)
    ax.set_ylabel('Entropy correlation btw\npatient and subtype B',
                  fontsize=fs)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.yaxis.set_tick_params(labelsize=fs)
    ax.grid(True)

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

