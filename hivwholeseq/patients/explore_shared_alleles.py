# vim: fdm=marker
'''
author:     Fabio Zanini
date:       02/07/14
content:    Explore parallel evolution in different patients by studying shared
            allele frequency trajectories.
'''
# Modules
import os
import argparse
from operator import itemgetter
import numpy as np
from Bio import SeqIO, AlignIO
import matplotlib.pyplot as plt



# Functions
def classify_sites(afts):
    '''Classify sites: 0 = conserved, 1 = polymorphic, 2 = sweeping'''
    site_class = np.zeros((afts.shape[0], afts.shape[1]), int)
    for i, aft in enumerate(afts):
        for j, aftpos in enumerate(aft):
            aftpos = np.vstack(aftpos)
            ia0 = aftpos[:, 0].argmax()
            if set((aftpos > 0.5).nonzero()[0]) - set([ia0]):
                site_class[i, j] = 2
    
            # FIXME: make this dependent on depth
            elif set((aftpos > 0.01).nonzero()[0]) - set([ia0]):
                site_class[i, j] = 1

    return site_class


def get_distance_class(site_class):
    '''Get distance matrix given a site class matrix'''
    d = np.zeros((site_class.shape[0], site_class.shape[0]), float)
    for i in xrange(d.shape[0]):
        for j in xrange(d.shape[0]):
            d[i, j] = ((site_class[i] != 0) != (site_class[j] != 0)).mean()
    return d


def reshape_single_aft(aft):
    '''Reshape a single aft (it's square now)'''
    aftn = np.ma.zeros((aft.shape[0], aft.shape[1], aft[0, 0].shape[0]), float)
    for i, af in enumerate(aft):
        aftn[i] = np.vstack(af)
    aftn.mask = aftn == -1
    return aftn



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get shared allele trajectories',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--regions', nargs='+', required=True,
                        help='Regions to analyze (e.g. F1 V3)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')

    args = parser.parse_args()
    regions = args.regions
    VERBOSE = args.verbose

    for region in regions:
        if VERBOSE >= 1:
            print region

        # Low-coverage regions are bytecoded by -1
        from hivwholeseq.patients.filenames import root_patient_folder
        data = np.load(root_patient_folder+'all/aft_shared_'+region+'.npz')
        afts = data['afts']
        pnames = data['pnames']
        depthmaxs = data['depthmaxs']
        n_patients = len(pnames)

        # Classes by sweeping, polymorphic, conserved
        site_class = classify_sites(afts)        
        d = get_distance_class(site_class)

        from scipy.cluster.hierarchy import linkage, dendrogram
        Z = linkage(d)
        
        # Plot as tree and image
        fig, axs = plt.subplots(1, 2, figsize=(16, 5), gridspec_kw={'width_ratios': [1, 4]})
        dg = dendrogram(Z, orientation='right', ax=axs[0], leaf_label_func=lambda x: '')
        axs[0].set_xticklabels('')
        ax = axs[1]
        ax.imshow(site_class[dg['leaves']][::-1], interpolation='nearest', aspect='auto')
        ax.set_yticks(np.arange(n_patients))
        ax.set_yticklabels(pnames[dg['leaves']][::-1].tolist())
        fig.suptitle(region)
        
        plt.tight_layout(rect=(0, 0, 1, 0.96))
        plt.ion()
        plt.show()
        
        ## Parallel sweeps in different patients        
        #is_sweeping = site_class == 2
        #h = np.bincount(is_sweeping.sum(axis=0))
        #fig, ax = plt.subplots()
        #plt.plot(np.arange(len(h)), h, lw=2, label='Data, '+region)
        #from scipy.stats import poisson
        #mu = curve_fit(poisson.pmf, np.arange(len(h)), 1.0 * h / h.sum(), p0=1)[0][0]
        #hth = poisson.pmf(np.arange(len(h)), mu)
        #plt.plot(np.arange(len(hth)), hth * h.sum(), lw=2, c='k',
        #         label='mu = '+'{:1.1e}'.format(mu)+', '+region)
        #
        #plt.grid()
        #plt.ylabel('# sweeping sites')
        #plt.xlabel('# patients')
        #plt.yscale('log')
        #plt.legend(loc=1)
        #
        #plt.ion()
        #plt.show()
