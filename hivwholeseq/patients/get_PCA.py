# vim: fdm=marker
'''
author:     Fabio Zanini
date:       01/12/14
content:    Perform PCA on the data to study population structure.
'''
# Modules
import sys
import os
import argparse
from itertools import izip
import numpy as np
from numpy import linalg
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.samples import SamplePat
from hivwholeseq.argparse_utils import RoiAction
import hivwholeseq.plot_utils




# Functions
def build_msa(htseqs, VERBOSE=0):
    '''Build multiple sequence alignment from cluster of haplotypes'''
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.Alphabet.IUPAC import ambiguous_dna
    
    seqs = [SeqRecord(Seq(seq, ambiguous_dna),
                      id='#'+str(i),
                      name='#'+str(i))
            for i, seq in enumerate(htseqs)]

    from hivwholeseq.mapping_utils import align_muscle
    ali = align_muscle(*seqs, sort=True)

    return ali


def check_main_splitting_sites(alibin, hft):
    '''Check which mutations are expected to make the difference in the PCA'''
    pass




# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Perform PCA on the data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--roi', required=True, action=RoiAction,
                        help='Region of interest (e.g. F1 300 350)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot local haplotype trajectories')

    args = parser.parse_args()
    pnames = args.patients
    roi = args.roi
    VERBOSE = args.verbose
    use_plot = args.plot

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]

    for pname, patient in patients.iterrows():

        patient = Patient(patient)
        patient.discard_nonsequenced_samples()

        if VERBOSE >= 1:
            print patient.name, roi
    
        if VERBOSE >= 2:
            print 'Get haplotype trajectories'
        try:
            (ht, indt, htseqs) = patient.get_region_count_trajectories(roi[0],
                                                                    VERBOSE=VERBOSE)
        except IOError:
            (ht, indt, htseqs) = patient.get_local_haplotype_count_trajectories(roi,
                                                                    VERBOSE=VERBOSE)
    
        if VERBOSE >= 2:
            print 'Align haplotypes and delete rare ones'
        indht = (ht > 5).any(axis=0)
        ht = ht[:, indht]
        htseqs = htseqs[indht]
        
        # Eliminate time points that are empty after this filter
        ind_keep = ht.any(axis=1)
        indt = indt[ind_keep]
        ht = ht[ind_keep]

        # Prepare data structures, both the matrix with sequences and the hft
        hft = (1.0 * ht.T / ht.sum(axis=1)).T
        alim = np.array(build_msa(htseqs))
        alibin = alim != alim[ht[0].argmax()]

        # NOTE: this function finds leads for the mutations that SHOULD matter,
        # as a sanity check for the later diagonalization
        check_main_splitting_sites(alibin, hft)
    
        if VERBOSE >= 2:
            print 'Restrict to polymorphic sites'
        # FIXME: use a better criterion
        ind_poly = (alibin.sum(axis=0) > 5)
        y = np.array(alibin[:, ind_poly], float)
    
        # Weight with the frequencies (from all time points??) and subtract mean
        # FIXME: decide on these weights!
        # for now normalize each time point by the coverange and weight equally
        w = (hft.sum(axis=0))**(0)
        yw = (w * y.T).T
        ym = yw - yw.mean(axis=0)
    
        if VERBOSE >= 2:
            print 'Diagonalize on covariance matrix'
        try:
            # NOTE: V[i] is the i-th eigenvector of the covariance matrix y.T x y, i.e.
            # V[i, j] is the projection of site j onto eigenvector i. For a sequence s,
            # its projection onto V[i] is sum_j (s[j] * V[i, j]), i.e. s V[i]. This is
            # also equal to U[i].
            U, s, V = linalg.svd(ym, full_matrices=False)

        # FIXME: there is an uncaught ValueError from DLASCL (parameter n4) in numpy
        except np.linalg.LinAlgError:
            if VERBOSE >= 1:
                print 'Diagonalization did not converge'
            continue
    
        if use_plot:
            # Divide the points by the weight of the haplotype, to fan them out a bit
            M3 = U[:, :3].T / w
            tis = ht.argmax(axis=0)
            colors = cm.jet(1.0 * ht.argmax(axis=0) / len(indt))
    
            # 2D plot
            fig, axs = plt.subplots(1, 2, figsize=(12, 6))
            ax = axs[0]
            ax.scatter(M3[0], M3[1], color=colors)
            ax.set_xlabel('PC1')
            ax.set_ylabel('PC2')
            ax.set_title(', '.join([patient.code] + map(str, roi)))
            ax.grid(True)
        
            # Haplotype trajectories
            from hivwholeseq.patients.get_local_haplotypes import plot_haplotype_frequencies
            plot_haplotype_frequencies(patient.times[indt], hft,
                                       figax=(fig, axs[1]),
                                       title=patient.name+', '+' '.join(map(str, roi)))
            #axs[1].set_yscale('logit')

            # 3D plot
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(M3[0], M3[1], M3[2], color=colors)
            ax.set_xlabel('PC1')
            ax.set_ylabel('PC2')
            ax.set_zlabel('PC3')
            ax.set_title(', '.join([patient.code] + map(str, roi)))


    if use_plot:
        plt.ion()
        plt.show()
    
