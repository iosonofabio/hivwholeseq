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

from hivwholeseq.generic_utils import mkdirs


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

