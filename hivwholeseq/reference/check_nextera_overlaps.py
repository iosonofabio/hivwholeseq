# vim: fdm=marker
'''
author:     Fabio Zanini
date:       15/12/13
content:    Check the bias in allele frequencies in overlaps in the Nextera
            dataset.
'''
# Modules
from itertools import izip
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.utils.miseq import alpha
from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.sequencing.filenames import get_allele_frequencies_filename
from hivwholeseq.sequencing.check_overlaps import get_overlapping_fragments, \
        get_overlap
from hivwholeseq.utils.mapping import align_muscle

import Tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Circle


# Globals
seq_run = 'Lina_nextera'
dataset = MiSeq_runs[seq_run]
data_folder = dataset['folder']
VERBOSE = 3
adaIDs = ['01', '02']



# Classes
class PointHighlighter(object):
    '''Functor to highlight points in both plots as a callback'''
    def __init__(self, data, axs, index_map=None):

        self.highlighted = []
        self.data = data
        self.axs = axs

        if index_map is not None:
            self.index_map = index_map
        else:
            if data[0][0].shape != data[1][0].shape:
                raise ValueError('The two datasets must be of the same size if there is no index_map')
            self.index_map = 2 * [range(data[0][0].shape[-1])]


    def __call__(self, event):
        if (not event.inaxes) or (event.inaxes not in self.axs):
            return
    
        # The first axis, data, index_map is always the clicked one
        if not self.axs.index(event.inaxes):
            axs = self.axs
            data = self.data
            index_map = self.index_map
        else:
            axs = self.axs[::-1]
            data = self.data[::-1]
            index_map = self.index_map[::-1]

        pos_click = (event.xdata, event.ydata)
        data_click = data[0]
        data_other = data[1]

        #             --------------X-------------        --------------Y------------
        ds = np.sqrt((data_click[0] - pos_click[0])**2 + (data_click[1] - pos_click[1])**2)
        ip_click_bas, ip_click_pos = np.unravel_index(np.argmin(ds), data_click[0].shape)
        point_click = (data_click[0][ip_click_bas, ip_click_pos],
                       data_click[1][ip_click_bas, ip_click_pos])

        ip_other_pos = index_map[0][ip_click_pos]
        if ip_other_pos is None:
            ip_other_bas = None
        else:
            ip_other_bas = ip_click_bas
            point_other = (data_other[0][ip_other_bas, ip_other_pos],
                           data_other[1][ip_other_bas, ip_other_pos])

        if self.highlighted:
            (ips_hi, points_hi, annos_hi) = self.highlighted
            points_hi[0].remove()
            points_hi[1].remove()
            annos_hi[0].remove()
            annos_hi[1].remove()
            self.highlighted = []              
            if ips_hi[0] == (ip_click_bas, ip_click_pos):
                plt.draw() 
                return

        ips_hi = []
        points_hi = []
        annos_hi = []
        axs[0].hold(True)
        scp = axs[0].scatter(point_click[0], point_click[1],
                             s=150, color='none', edgecolor='r', lw=2)
        anno = axs[0].annotate(str(ip_click_pos)+alpha[ip_click_bas],
                               (0.05, 0.9), xycoords="axes fraction",
                               bbox=dict(fc='w', boxstyle='round'))
        annos_hi.append(anno)
        points_hi.append(scp)
        ips_hi.append((ip_click_bas, ip_click_pos))
        axs[0].hold(False)

        axs[1].hold(True)
        if ip_other_pos is not None:
            scp = axs[1].scatter(point_other[0], point_other[1],
                                 s=150, color='none', edgecolor='r', lw=2)
            anno = axs[1].annotate(str(ip_other_pos)+alpha[ip_other_bas],
                                   (0.05, 0.9), xycoords="axes fraction",
                                   bbox=dict(fc='w', boxstyle='round'))
        else:
            scp = axs[1].text(0.5, 0.5, 'X', color='r', fontsize=80, transform=axs[1].transAxes)
            anno = axs[1].annotate('Not found',
                                   (0.05, 0.9), xycoords="axes fraction",
                                   bbox=dict(fc='w', boxstyle='round'))
        annos_hi.append(anno)
        ips_hi.append((ip_other_bas, ip_other_pos))
        points_hi.append(scp)
        axs[1].hold(False)

        self.highlighted = (ips_hi, points_hi, annos_hi)

        plt.draw() 
        return



# Functions
def check_overlap_allele_frequencies(data_folder, adaID, frag1, frag2, overlap,
                                     VERBOSE=0, ax=None):
    '''Check biases in allele frequencies in the overlap'''
    (start_s2, end_s1, ali) = overlap

    # Reduce the allele counts (fwd and rev have different status on most overlaps,
    # because of the uneven coverage)
    nu1 = np.load(get_allele_frequencies_filename(data_folder, adaID, frag1))[:,start_s2:]
    nu2 = np.load(get_allele_frequencies_filename(data_folder, adaID, frag2))[:,:end_s1] 

    # FIXME
    if nu1.shape != nu2.shape:
        return

    print nu1.shape

    # Print table of called polymorphisms
    thre = 3e-3
    print 'adaID', adaID, frag1, frag2, 'polymorphism matrix (NO | YES), threshold = '+\
            '{:1.1e}'.format(thre)
    print 3 * ' ', '|', '{:^10s}'.format(frag1)
    print 15 * '-'
    print 3 * ' ', '|', \
            '{:3d}'.format(((nu1 < 1e-6) & (nu2 < 1e-6)).sum()), '|', \
            '{:3d}'.format(((nu1 > thre) & (nu2 < 1e-6)).sum())
    print '{:3s}'.format(frag2), '+'+(5*'-')+'+'+(4*'-')
    print 3 * ' ', '|', \
            '{:3d}'.format(((nu1 < 1e-6) & (nu2 > thre)).sum()), '|', \
            '{:3d}'.format(((nu1 > thre) & (nu2 > thre)).sum())
         

    # Plot scatter
    xdata = np.abs(nu1 - 1e-5)
    ydata = np.abs(nu2 - 1e-5)
    sc = ax.scatter(xdata, ydata, s=30,
                    c=cm.jet([int(255.0 * i / len(nu1.T)) for i in xrange(len(nu1.T))]))

    # Plot diagonal
    ax.plot([1e-7, 2], [1e-7, 2], lw=1, c='k', ls='--')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.7e-5, 1.2)
    ax.set_ylim(0.7e-5, 1.2)

    return (xdata, ydata)


def load_data():
    overlaps = {}
    pairs = get_overlapping_fragments(fragments)
    for adaID in adaIDs:
        overlaps[adaID] = {}
        for (frag1, frag2) in pairs:
            overlap = get_overlap(data_folder, adaID, frag1, frag2, VERBOSE=VERBOSE)
            overlaps[adaID][(frag1, frag2)] = overlap
    return overlaps


def make_index_map(overlaps):
    overlap1 = overlaps[0]
    overlap2 = overlaps[1]

    ali1 = overlap1[2]
    ali2 = overlap2[2]

    ali11 = align_muscle(ali1[0], ali2[0])

    index_map = [[], []]
    pos1 = 0
    pos2 = 0
    for a1, a2 in izip(*ali11):
        if (a1 != '-') and (a2 != '-'):
            index_map[0].append(pos2)
            index_map[1].append(pos1)
            pos1 += 1
            pos2 += 1
        elif a1 == '-':
            index_map[1].append(None)
            pos2 += 1
        else:
            index_map[0].append(None)
            pos1 += 1

    return index_map



# Script
if __name__ == '__main__':

    # FIXME: do for all
    fragments = ['F'+str(i) for i in xrange(1, 3)]

    frag1 = fragments[0] 
    frag2 = fragments[1] 

    fig, axs = plt.subplots(2, 1, figsize=(10, 14))
    fig.suptitle('Nextera dataset, '+frag1+'-'+frag2+' overlap') 
    axs[0].set_title('adaID 01')
    axs[1].set_title('adaID 02')
    axs[1].set_xlabel(r'$\nu_{'+frag1+r'}$', fontsize=20)
    axs[0].set_ylabel(r'$\nu_{'+frag2+r'}$', fontsize=20)
    axs[1].set_ylabel(r'$\nu_{'+frag2+r'}$', fontsize=20)

    overlap1 = get_overlap(data_folder, '01', frag1, frag2, VERBOSE=VERBOSE)
    (start_s2, end_s1, ali1) = overlap1
    nu11 = np.load(get_allele_frequencies_filename(data_folder, '01', frag1))[:,start_s2:]
    nu12 = np.load(get_allele_frequencies_filename(data_folder, '01', frag2))[:,:end_s1] 

    data1 = check_overlap_allele_frequencies(data_folder, '01', frag1, frag2, overlap1,
                                             VERBOSE=VERBOSE, ax=axs[0])

    overlap2 = get_overlap(data_folder, '02', frag1, frag2, VERBOSE=VERBOSE)
    (start_s2, end_s1, ali2) = overlap2
    nu21 = np.load(get_allele_frequencies_filename(data_folder, '02', frag1))[:,start_s2:][:, 4:]
    nu22 = np.load(get_allele_frequencies_filename(data_folder, '02', frag2))[:,:end_s1][:, 4:]

    data2 = check_overlap_allele_frequencies(data_folder, '02', frag1, frag2, overlap2,
                                             VERBOSE=VERBOSE, ax=axs[1])
    data2 = (data2[0], data2[1])

    index_map = make_index_map((overlap1, overlap2))

    highlighter = PointHighlighter((data1, data2), list(axs), index_map)
    plt.connect('button_press_event', highlighter)

    plt.tight_layout(rect=(0, 0, 1, 0.95))
    plt.ion()
    plt.show()
