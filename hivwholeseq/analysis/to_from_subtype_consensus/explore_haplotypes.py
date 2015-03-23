# vim: fdm=indent
'''
author:     Fabio Zanini
date:       18/03/15
content:    What part of the sequence space (tree, cluster) do haplotypes explore
            within a patient? This is relevant because if the patient were to infect
            somebody else at that point, those are the potential leaves for the
            updated tree.
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
from Bio.Seq import translate

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
from hivwholeseq.utils.sequence import get_coordinates_genomic_region

from hivwholeseq.analysis.to_from_subtype_consensus.random_haplotypes import (
    align_to_reference)



# Globals
pnames = ['20097', '15319', '9669', '15363', '15823', '15241', '15376']
regions = ['p17']




# Functions
def collect_data(pnames, regions, VERBOSE=0, plot=False):
    '''Collect data to study allele freqs around sweeps'''
    from hivwholeseq.reference import load_custom_reference

    regdata = []
    data = []

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    refseq = load_custom_reference('HXB2', 'gb')

    for region in regions:
        if VERBOSE >= 1:
            print region

        if VERBOSE >= 2:
            print 'Get subtype alignment'
        ali = get_subtype_reference_alignment(region, VERBOSE=VERBOSE)
        alim = np.array(ali)

        if VERBOSE >= 2:
            print 'Get subtype consensus (for checks only)'
        consrec = get_subtype_reference_alignment_consensus(region, VERBOSE=VERBOSE)
        conss = ''.join(consrec)

        if VERBOSE >= 2:
            print 'Get reference sequence'
        location = get_coordinates_genomic_region(refseq, region)
        start = location.nofuzzy_start
        end = location.nofuzzy_end

        regdata.append({'name': region,
                        'ali': ali, 'alim': alim,
                        'consensus': consrec,
                        'location': (start, end),
                        'start': start, 'end': end,
                        'L': end - start})

        for pname, patient in patients.iterrows():
            if VERBOSE >= 1:
                print pname, region
            patient = Patient(patient)

            if VERBOSE >= 2:
                print 'Get haplotype count'
            hct, ind, alim = patient.get_haplotype_count_trajectory(region)
            alism = np.array(map(''.join, alim), 'S'+str(len(alim[0])))

            # Exclude time points without counts
            ii = (hct > 0).any(axis=1)
            if not ii.sum():
                continue

            hct = hct[ii]
            ind = ind[ii]
            times = patient.times[ind]

            if VERBOSE >= 2:
                print 'Align to reference'
            alisma = align_to_reference(conss, alism, VERBOSE=VERBOSE)
            alima = np.array([np.fromstring(seq, 'S1') for seq in alisma])

            if VERBOSE >= 2:
                print 'Get frequencies'
            hft = (1.0 * hct.T / hct.sum(axis=1)).T

            if VERBOSE >= 2:
                print 'Get initial consensus'
            consm0 = alima[hft[0].argmax()]

            data.append({'region': region, 'pcode': patient.code,
                         'hft': hft, 'times': times, 'ali': alima})

        data = pd.DataFrame(data)
        regdata = pd.DataFrame(regdata)
        regdata.set_index('name', drop=False, inplace=True)

        return {'data': data, 'regdata': regdata}


def get_pairwise_distance_alignments(alim1, alim2, for_cycle=False):
    '''Get pairwise distance between two sequence matrices
    
    Parameters:
       for_cycle (bool): use a python iterator instead of BIG numpy tensors
    '''
    if alim1.shape[-1] != alim2.shape[-1]:
        raise ValueError('The two alignments must have the same length')

    N1 = alim1.shape[0]
    N2 = alim2.shape[0]
    L = alim1.shape[-1]

    if not for_cycle:
        a = np.zeros((N1, N2, L), 'S1')
        a[:] = alim2
    
        b = np.zeros((N2, N1, L), 'S1')
        b[:] = alim1

        d = (a != b.swapaxes(0, 1)).sum(axis=-1)
        return d

    # Iterate over the shorter axis
    d = np.zeros((N1, N2), int)
    if N1 <= N2:
        for i, seq1 in enumerate(alim1):
            d[i] = (seq1 != alim2).sum(axis=-1)
    else:
        for i, seq2 in enumerate(alim2):
            d[:, i] = (seq2 != alim1).sum(axis=-1)

    return d


def take_close_subalignment(alim, alipam, percentile=0.5, VERBOSE=0):
    '''Take a subalignment of the subtype that is close to the patients'''
    # Distance matrix
    d = get_pairwise_distance_alignments(alipam, alim)

    # Take the closest seqs to any patient haplotype
    q = np.percentile(d, percentile, axis=-1)
    ind = np.unique((d.T < q).T.nonzero()[1])

    return alim[ind], ind


def plot_clusterforce(v, dcon):
    '''Plot force field clustering'''
    from matplotlib import cm
    import matplotlib.pyplot as plt

    # Plot the force field and the scatter
    fig, ax = plt.subplots()

    colors = cm.jet(1.0 * dcon / dcon.max())
    ax.scatter([0], [0], s=200, edgecolor='k', facecolor='none', lw=2, zorder=-1)
    ax.scatter(v[:, 0], v[:, 1], s=40, c=colors)
    ax.grid(True)
    ax.set_xlim(-1.04*np.abs(v[:, 0]).max(), 1.04*np.abs(v[:, 0]).max())
    ax.set_ylim(-1.04*np.abs(v[:, 1]).max(), 1.04*np.abs(v[:, 1]).max())
    sm = plt.cm.ScalarMappable(cmap=cm.jet, norm=plt.Normalize(vmin=0, vmax=dcon.max()))
    sm.set_array(dcon)
    cb = plt.colorbar(sm)
    cb.set_label('Hamming distance from consensus', rotation=270, labelpad=40)
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    plt.tight_layout()

    return fig, ax


def draw_tree(tree, label_func=str, do_show=True, show_confidence=True,
         # For power users
         fontsize=None, linewidth=None,
         axes=None, branch_labels=None, *args, **kwargs):
    """Plot the given tree using matplotlib (or pylab).

    The graphic is a rooted tree.

    Additional keyword arguments passed into this function are used as pyplot
    options. The input format should be in the form of:
    pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict), or
    pyplot_option_name=(dict).

    Example using the pyplot options 'axhspan' and 'axvline':

    >>> Phylo.draw(tree, axhspan=((0.25, 7.75), {'facecolor':'0.5'}),
    ...     axvline={'x':'0', 'ymin':'0', 'ymax':'1'})

    Visual aspects of the plot can also be modified using pyplot's own functions
    and objects (via pylab or matplotlib). In particular, the pyplot.rcParams
    object can be used to scale the font size (rcParams["font.size"]) and line
    width (rcParams["lines.linewidth"]).

    :Parameters:
        label_func : callable
            A function to extract a label from a node. By default this is str(),
            but you can use a different function to select another string
            associated with each node. If this function returns None for a node,
            no label will be shown for that node.
        do_show : bool
            Whether to show() the plot automatically.
        show_confidence : bool
            Whether to display confidence values, if present on the tree.
        axes : matplotlib/pylab axes
            If a valid matplotlib.axes.Axes instance, the phylogram is plotted
            in that Axes. By default (None), a new figure is created.
        branch_labels : dict or callable
            A mapping of each clade to the label that will be shown along the
            branch leading to it. By default this is the confidence value(s) of
            the clade, taken from the ``confidence`` attribute, and can be
            easily toggled off with this function's ``show_confidence`` option.
            But if you would like to alter the formatting of confidence values,
            or label the branches with something other than confidence, then use
            this option.
    """

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        try:
            import pylab as plt
        except ImportError:
            from Bio import MissingPythonDependencyError
            raise MissingPythonDependencyError(
                "Install matplotlib or pylab if you want to use draw.")

    import matplotlib.collections as mpcollections

    # Arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []

    # Options for displaying branch labels / confidence
    def conf2str(conf):
        if int(conf) == conf:
            return str(int(conf))
        return str(conf)
    if not branch_labels:
        if show_confidence:
            def format_branch_label(clade):
                if hasattr(clade, 'confidences'):
                    # phyloXML supports multiple confidences
                    return '/'.join(conf2str(cnf.value)
                                    for cnf in clade.confidences)
                if clade.confidence:
                    return conf2str(clade.confidence)
                return None
        else:
            def format_branch_label(clade):
                return None
    elif isinstance(branch_labels, dict):
        def format_branch_label(clade):
            return branch_labels.get(clade)
    else:
        assert callable(branch_labels), \
            "branch_labels must be either a dict or a callable (function)"
        format_branch_label = branch_labels

    # Options for line width, font size, etc.
    if fontsize is None:
        fontsize = plt.rcParams['font.size']
    if linewidth is None:
        linewidth = plt.rcParams['lines.linewidth']


    # Layout
    def get_x_positions(tree):
        """Create a mapping of each clade to its horizontal position.

        Dict of {clade: x-coord}
        """
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
        return depths

    def get_y_positions(tree):
        """Create a mapping of each clade to its vertical position.

        Dict of {clade: y-coord}.
        Coordinates are negative, and integers for tips.
        """
        maxheight = tree.count_terminals()
        # Rows are defined by the tips
        heights = dict((tip, maxheight - i)
                       for i, tip in enumerate(reversed(tree.get_terminals())))

        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = (heights[clade.clades[0]] +
                              heights[clade.clades[-1]]) / 2.0

        if tree.root.clades:
            calc_row(tree.root)
        return heights

    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)
    # The function draw_clade closes over the axes object
    if axes is None:
        fig, axes = plt.subplots()
    elif not isinstance(axes, plt.matplotlib.axes.Axes):
        raise ValueError("Invalid argument for axes: %s" % axes)

    def draw_clade_lines(use_linecollection=False, orientation='horizontal',
                         y_here=0, x_start=0, x_here=0, y_bot=0, y_top=0,
                         color='black', lw='.1'):
        """Create a line with or without a line collection object.

        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """
        if (use_linecollection is False and orientation == 'horizontal'):
            axes.hlines(y_here, x_start, x_here, color=color, lw=lw)
        elif (use_linecollection is True and orientation == 'horizontal'):
            horizontal_linecollections.append(mpcollections.LineCollection(
                [[(x_start, y_here), (x_here, y_here)]], color=color, lw=lw),)
        elif (use_linecollection is False and orientation == 'vertical'):
            axes.vlines(x_here, y_bot, y_top, color=color)
        elif (use_linecollection is True and orientation == 'vertical'):
            vertical_linecollections.append(mpcollections.LineCollection(
                [[(x_here, y_bot), (x_here, y_top)]], color=color, lw=lw),)

    def draw_clade(clade, x_start, color, lw=1, fontsize='small'):
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        # phyloXML-only graphics annotations
        if hasattr(clade, 'color') and clade.color is not None:
            color = clade.color.to_hex()
        if hasattr(clade, 'width') and clade.width is not None:
            lw = clade.width * plt.rcParams['lines.linewidth']
        # Draw a horizontal line from start to here
        draw_clade_lines(use_linecollection=True, orientation='horizontal',
                         y_here=y_here, x_start=x_start, x_here=x_here,
                         color=color, lw=lw)
        # Add node/taxon labels
        label = label_func(clade)
        if label not in (None, clade.__class__.__name__):
            axes.text(x_here, y_here, ' %s' %
                      label, verticalalignment='center')
        # Add label above the branch (optional)
        conf_label = format_branch_label(clade)
        if conf_label:
            axes.text(0.5 * (x_start + x_here), y_here, conf_label,
                      fontsize=fontsize, horizontalalignment='center')
        if clade.clades:
            # Draw a vertical line connecting all children
            y_top = y_posns[clade.clades[0]]
            y_bot = y_posns[clade.clades[-1]]
            # Only apply widths to horizontal lines, like Archaeopteryx
            draw_clade_lines(use_linecollection=True, orientation='vertical',
                             x_here=x_here, y_bot=y_bot, y_top=y_top,
                             color=color, lw=lw)
            # Draw descendents
            for child in clade:
                draw_clade(child, x_here, color, lw)

    draw_clade(tree.root, 0, 'k', lw=linewidth, fontsize=fontsize)

    # If line collections were used to create clade lines, here they are added
    # to the pyplot plot.
    for i in horizontal_linecollections:
        axes.add_collection(i)
    for i in vertical_linecollections:
        axes.add_collection(i)

    # Aesthetics

    if hasattr(tree, 'name') and tree.name:
        axes.set_title(tree.name, fontsize=fontsize)
    axes.set_xlabel('branch length', fontsize=fontsize)
    axes.set_ylabel('taxa', fontsize=fontsize)
    # Add margins around the tree to prevent overlapping the axes
    xmax = max(x_posns.values())
    axes.set_xlim(-0.05 * xmax, 1.25 * xmax)
    # Also invert the y-axis (origin at the top)
    # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
    axes.set_ylim(max(y_posns.values()) + 0.8, 0.2)

    # Parse and process key word arguments as pyplot options
    for key, value in kwargs.items():
        try:
            # Check that the pyplot option input is iterable, as required
            [i for i in value]
        except TypeError:
            raise ValueError('Keyword argument "%s=%s" is not in the format '
                             'pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict),'
                             ' or pyplot_option_name=(dict) '
                             % (key, value))
        if isinstance(value, dict):
            getattr(plt, str(key))(**dict(value))
        elif not (isinstance(value[0], tuple)):
            getattr(plt, str(key))(*value)
        elif (isinstance(value[0], tuple)):
            getattr(plt, str(key))(*value[0], **dict(value[1]))

    if do_show:
        plt.show()

    return {'positions': {'x': x_posns, 'y': y_posns},
            'ax': axes,
           }




# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Mutations away from/towards subtype',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+', default=pnames,
                        help='Patients to analyze')
    parser.add_argument('--regions', nargs='+', default=regions,
                        help='Genomic regions (e.g. V3 IN)')
    parser.add_argument('--verbose', type=int, default=2,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot results')

    args = parser.parse_args()
    pnames = args.patients
    regions = args.regions
    VERBOSE = args.verbose
    use_plot = args.plot


    res = collect_data(pnames, regions, VERBOSE=VERBOSE, plot=use_plot)
    data = res['data']
    regdata = res['regdata']

   
    if VERBOSE >= 1:
        'Analyze data'
    for region, datum in data.groupby('region'):
        if VERBOSE >= 1:
            print region


        ali = regdata.loc[region, 'ali']
        alim = regdata.loc[region, 'alim']
        consm = np.array(regdata.loc[region, 'consensus'], 'S1')
        
        alipam = np.concatenate(datum['ali'])

        if VERBOSE >= 2:
            print 'Take subalignment'
        alisubm, indsub = take_close_subalignment(alim, alipam)

        # Visualize subtype alignment in various ways
        # Force field clustering
        from clusterforce.clustering import cluster_force, position_sequence
        res = cluster_force(alisubm, method='BFGS-jac', plot=False)
        v = res['x']

        dcon = (alisubm != consm).sum(axis=-1)

        # For each patient, add the points according to the frozen force field
        for pcode, datumpa in datum.groupby('pcode'):
            fig, ax = plot_clusterforce(v, dcon)

            alipam = datumpa.iloc[0]['ali']
            hft = datumpa.iloc[0]['hft']
            times = datumpa.iloc[0]['times']

            us = []
            for seqpa, hf in izip(alipam, hft.T):
                u = position_sequence(seqpa, v, alisubm, e1e2=res['e1e2'])

                # Take the time of maximal frequency
                it = hf.argmax()
                time = times[it]
                nu = hf[it]

                numin = 0.01
                if nu < numin:
                    continue

                label = '{:.1G}, {:d} m'.format(nu, int(time / 30.5))
                color = cm.jet(1.0 * time / times.max())
                color = tuple(list(color[:-1]) + [0.4])

                size = 10 + 290 * (np.log(nu) - np.log(numin)) / (0.0 - np.log(numin))

                ax.scatter(u[0], u[1],
                           s=size,
                           marker='s',
                           edgecolor='k', facecolor=color, lw=1)
                #ax.text(u[0], u[1], label, fontsize=8)

                us.append(u)

            ax.set_title(pcode)

            uall = np.vstack([v, us])
            xmin = uall[:, 0].min()
            xmax = uall[:, 0].max()
            xspan = xmax - xmin
            ymin = uall[:, 1].min()
            ymax = uall[:, 1].max()
            yspan = ymax - ymin
            ax.set_xlim(xmin - 0.04 * xspan, xmax + 0.04 * xspan)
            ax.set_ylim(ymin - 0.04 * yspan, ymax + 0.04 * yspan)


        # Phylogenetic tree with ancestral sequences
        from hivwholeseq.cross_sectional.get_subtype_reference_alignment_tree import (
            get_subtype_reference_alignment_tree)
        tree = get_subtype_reference_alignment_tree(region, format='json')
        # Prune to subalignment
        subnames = set(ali[i].name for i in indsub)
        leaves = tree.get_terminals()
        for leaf in leaves:
            if leaf.name not in subnames:
                tree.prune(leaf)
            else:
                subnames.remove(leaf.name)

        # Remove labels from internal nodes
        for node in tree.get_nonterminals():
            node.name = None

        # Add distance from consensus
        for node in tree.get_terminals() + tree.get_nonterminals():
            if not hasattr(node, 'sequence'):
                node.color = [0, 0, 0]
            else:
                seq = np.fromstring(node.sequence, 'S1')
                # FIXME: this requires alignment!
                dseq = (consm != seq).sum()
                node.distance = dseq
                node.color = [int(255.0 * c) for c in cm.jet(1.0 * dseq / dcon.max())[:3]]

        # Get sequences and node pointers
        nodes, seqtree = zip(*[(node, np.fromstring(node.sequence, 'S1'))
                               for node in tree.get_terminals() + tree.get_nonterminals()
                               if node.sequence is not None])
        seqtree = np.array(seqtree)

        # For each patient, add the points
        for pcode, datumpa in datum.groupby('pcode'):
            if pcode != 'p8':
                continue
            # Plot the tree
            fig, ax = plt.subplots()
            datap = draw_tree(tree, axes=ax, show_confidence=False, fontsize=8,
                              linewidth=1.5, do_show=False,
                              label_func=lambda x: None,
                             )
            nodespos = datap['positions']

            alipam = datumpa.iloc[0]['ali']
            hft = datumpa.iloc[0]['hft']

            d = get_pairwise_distance_alignments(alipam, seqtree)
            imin = d.argmin(axis=1)
            iminu = np.unique(imin)
            labels = []
            for i in iminu:
                hfs = hft[:, imin == i]
                label = '\n'.join([' '.join(['{:.1G}'.format(nu) for nu in hf])
                                   for hf in hfs.T])
                labels.append(label)

            for i in xrange(len(iminu)):
                node = nodes[iminu[i]]
                x = nodespos['x'][node]
                y = nodespos['y'][node]
                ax.scatter(x, y, s=100, marker='s', edgecolor='k', facecolor='none',
                           lw=3, zorder=100)

            ax.set_title(pcode)
            ax.grid(True)



        plt.ion()
        plt.show()
