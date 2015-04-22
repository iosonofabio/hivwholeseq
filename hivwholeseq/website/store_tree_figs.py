# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/12/14
content:    Plot the quality along the read pair.
'''
# Modules
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patient
from hivwholeseq.utils.sequence import align_muscle
from hivwholeseq.utils.tree import (build_tree_fasttree, filter_rare_leaves,
                                    correct_minimal_branches)

from hivwholeseq.website.filenames import get_tree_figure_filename
from hivwholeseq.paper_figures.plots import HIVEVO_colormap

# Globals
pnames = ['20097', '15363', '15823', '15376', '9669', '15107', '15241', '15034', '15319']
region = 'V3'
cutoff = 0.04

web_colormap = HIVEVO_colormap()



# Functions
def compress_data(tree, pname, region):
    '''Compress data for plots, discarding useless info'''
    data = []
    datum = {'tree': tree,
             'pname': pname,
             'region': region}
    data.append(datum)

    return data


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
        tree.root.branch_length = 0.001

        times = sorted(set([leaf.DSI for leaf in tree.get_terminals()]))

        assign_color(tree, attrname=color_by, cmap = web_colormap)
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
        ax.set_xlabel('Genetic distance [changes / site]', fontsize=16, labelpad=10)

        # Add circles to the leaves
        (x, y, s, c,cs) = zip(*data_circles)
        ax.scatter(x, y, s=s, c=c, edgecolor=cs, zorder=2)
        ax.set_xlim(-0.04 * maxdepth, 1.04 * maxdepth)

        # Draw a "legend" for sizes
        datal = [{'hf': 0.05, 'label': '5%'},
                 {'hf': 0.20, 'label': '20%'},
                 {'hf': 1.00, 'label': '100%'}]
        ax.text(0.98 * maxdepth, 0.01 * ax.get_ylim()[0],
                'Haplotype frequency:', fontsize=16, ha='right')
        for idl, datuml in enumerate(datal):
            r = rfun(datuml['hf'])
            y = (0.07 + 0.07 * idl) * ax.get_ylim()[0]
            ax.scatter(0.85 * maxdepth, y, s=r,
                       facecolor='k',
                       edgecolor='none')
            ax.text(0.98 * maxdepth, y + 0.02 * ax.get_ylim()[0],
                    datuml['label'], ha='right',
                    fontsize=14)

        # Draw legend for times
        datal = [{'time': t,'color': [tone for tone in web_colormap(1.0 * t / max(times))[:3]],
        'colorstroke': [tone * 0.7 for tone in web_colormap(1.0 * t / max(times))[:3]]} for t in times]
        xtext = 0
        ytext = (0.93 - 0.06 * min(8, len(datal))) * ax.get_ylim()[0]
        ax.text(0.01 * maxdepth, ytext, 'Time:', fontsize=16)
        for ico, datuml in enumerate(datal):
            ytext += 0.06 * (ax.get_ylim()[0] - ax.get_ylim()[1])

            if ico == 8:
                ytext -= 4 * 0.06 * (ax.get_ylim()[0] - ax.get_ylim()[1])
                xtext += 0.28 * maxdepth

            ax.scatter(xtext, ytext, s=rfun(0.5),
                       facecolor=datuml['color'],
                       edgecolor=datuml['colorstroke'])
            ax.text(xtext + 0.21 * maxdepth, ytext + 0.02 * ax.get_ylim()[0],
                    str(int(datuml['time'] / 30.5))+' months',
                    ha='right',
                    fontsize=14,
                   )

        # Draw scale bar
        xbar = (0.3 + 0.3 * (len(datal) >= 9)) * maxdepth
        ybar = 0.90 * ax.get_ylim()[0]
        lbar = 0.05 * maxdepth
        lbar_label = '{:.1G}'.format(lbar)
        lbar = float(lbar_label)
        ax.plot([xbar, xbar + lbar], [ybar, ybar], lw=4, c='k')
        ax.text(xbar + 0.5 * lbar, ybar + 0.08 * ax.get_ylim()[0],
                lbar_label, fontsize=14,
                ha='center')


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




# Script
if __name__ == '__main__':

    VERBOSE = 2

    for pname in pnames:
        print pname

        patient = load_patient(pname)
        tree = patient.get_local_tree(region)

        filter_rare_leaves(tree, cutoff)
        
        correct_minimal_branches(tree)

        data = compress_data(tree, patient.code, region)
            
        filename = get_tree_figure_filename(patient.code, region, format='svg')

        plot_haplotype_tree_example(data,
                                    VERBOSE=VERBOSE,
                                    savefig=filename)

        #plot_haplotype_tree_example(data,
        #                            VERBOSE=VERBOSE,
        #                            )
        #break
