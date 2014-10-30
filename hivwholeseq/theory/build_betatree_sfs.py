#!/usr/bin/env python
# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/10/14
content:    Precompute a few site frequency spectra from betatrees at different
            N and alpha, so we have them for HIV comparisons.
'''
# Modules
import os
import argparse
import numpy as np

from hivwholeseq.fork_cluster import fork_calculate_beta_SFS as fork_self



# Functions
def calculate_beta_SFS(VERBOSE=0, alpha=1, bins=20, N=100, save=True, ntrees=1000):
    '''Load sfs of direct bsc simulations'''
    import numpy as np
    from hivwholeseq.theory.betatree.src.sfs import SFS
    from hivwholeseq.theory.filenames import get_sfs_betatree_filename
    
    sfs_beta = SFS(sample_size=N, alpha=alpha)

    if VERBOSE >= 2:
        print 'Generating beta coalescent SFS with N = '+str(N)+' and alpha = '+str(alpha)
    sfs_beta.getSFS(ntrees=ntrees)

    sfs_beta.binSFS(mode='logit', bins=bins)
    if np.iterable(bins):
        sfs_beta.bin_center = np.sqrt(bins[1:] * bins[:-1])

    x = sfs_beta.bin_center
    y = sfs_beta.binned_sfs

    data = np.vstack([x, y]).T

    if save:
        fn_out = get_sfs_betatree_filename(N, alpha)
        np.savetxt(fn_out, data, header='# bin center\tSFS (density)', delimiter='\t')

    return data.T



# Script
if __name__ == '__main__':

    alphas = [1, 1.25, 1.5, 1.75]
    Ns = [300, 1000, 3000, 10000]

    parser = argparse.ArgumentParser(description='Get site frequency spectra',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--Ns', nargs='+', type=int, default=Ns,
                        help='Sample sizes')
    parser.add_argument('--alphas', nargs='+', type=float, default=alphas,
                        help='Sample sizes')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    Ns = args.Ns
    alphas = args.alphas
    submit = args.submit
    VERBOSE = args.verbose

    # Logit bins
    tbins = np.linspace(-7,7,21)
    bins = np.exp(tbins)/(1+np.exp(tbins))

    sfss = {}
    for N in Ns:
        for alpha in alphas:
            if submit:
                fork_self(alpha=alpha, N=N, VERBOSE=VERBOSE)
                continue

            sfss[(N, alpha)] = calculate_beta_SFS(VERBOSE=VERBOSE, alpha=alpha, N=N)
            
