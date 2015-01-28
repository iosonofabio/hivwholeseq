# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/01/15
content:    Correlate entropy in a subtype alignment with the fitness mutation
            profile from Al-Mawsawi et al 2014.
'''
# Modules
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha
from hivwholeseq.reference import load_custom_reference
from hivwholeseq.utils.sequence import get_coordinates_genomic_region
from hivwholeseq.patients.explore_patconserved_subtype import get_subtype_reference_alignment, get_ali_entropy
from hivwholeseq.reference_samples.mutational_profiling_NL43 import load_table


# Functions
def get_consensus_nuc(ali, positions):
    '''Get consensus nucleotide for each of those positions'''
    from collections import Counter

    cons = []
    for pos in positions:
        col = ali[:, pos]
        cons.append(Counter(col).most_common(1)[0][0])
    return cons



def take_minimal_cost(table_reg):
    '''Take only the minimal cost for positions that have several data points'''
    ind = np.zeros(table_reg.shape[0], bool)

    # NOTE: this algorithm is not optimized at all, but it should be fine
    pos_unique = np.unique(table_reg['Pos HXB2'])
    for pos in pos_unique:
        ind_pos = (table_reg['Pos HXB2'] == pos).nonzero()[0]
        if len(ind_pos) == 1:
            ind[ind_pos[0]] = True
            continue

        table_sub = table_reg.iloc[ind_pos]
        # NOTE: here numpy and pandas work differently
        i_maxF = np.array(table_sub['RC Index']).argmax()
        ind[ind_pos[i_maxF]] = True

    return table_reg.loc[ind]



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Correlate subtype entropy and fitness',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Regions to analyze (e.g. F1 p17)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot results')

    args = parser.parse_args()
    regions = args.regions
    VERBOSE = args.verbose
    plot = args.plot
 
    table = load_table(align_to_reference='HXB2')
    ref = load_custom_reference('HXB2', format='gb')

    for region in regions:
        ali = get_subtype_reference_alignment(region, VERBOSE=VERBOSE, refname='HXB2')

        posregion = get_coordinates_genomic_region(ref, region)
        pos_start = posregion.nofuzzy_start
        # NOTE: Pavel sometimes took the first subdomain (e.g. RT)
        pos_end = pos_start + ali.get_alignment_length()
        ind_reg = ((table['Pos HXB2'] >= pos_start) &
                   (table['Pos HXB2'] < pos_end))
    
        table_reg = table.loc[ind_reg.nonzero()[0]].copy()
        posali = table_reg['Pos HXB2'] - pos_start

        table_reg['Nuc HXB2'] = np.array(ref)[table_reg['Pos HXB2']]
        table_reg['Nuc subtype B'] = get_consensus_nuc(ali, posali)

        table_reg = take_minimal_cost(table_reg)
        posali = table_reg['Pos HXB2'] - pos_start

        # RC index is: frequency after 2 passages / initial frequency, so
        # the log should be fine for a fitness
        fitness = np.log10(table_reg['RC Index'])
        entropy = get_ali_entropy(ali, posali)

        from scipy.stats import spearmanr, pearsonr
        print 'Spearman R:', spearmanr(fitness, entropy)
        print 'Pearson r:', pearsonr(fitness, entropy)

        if plot:
            fig, ax = plt.subplots()
            ax.scatter(entropy, fitness)
            ax.set_xlabel('Entropy in subtype B')
            ax.set_ylabel('Fitness_10')
            ax.set_xscale('log')
            ax.set_xlim(1e-3, 2)
            ax.grid(True)
            ax.set_title(region)

            plt.tight_layout()
            plt.ion()
            plt.show()
