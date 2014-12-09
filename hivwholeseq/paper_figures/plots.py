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


def plot_rna_amplification_bias(data, title='', VERBOSE=0, savefig=False):
    '''Plot amplification bias of the RNA mixes'''
    #TODO
    pass
