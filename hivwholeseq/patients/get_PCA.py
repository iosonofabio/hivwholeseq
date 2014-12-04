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
from hivwholeseq.patients.patients import load_patient
from hivwholeseq.patients.samples import SamplePat
from hivwholeseq.argparse_utils import RoiAction
from hivwholeseq.patients.get_roi import get_fragmented_roi



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


# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Perform PCA on the data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--roi', required=True, action=RoiAction,
                        help='Region of interest (e.g. F1 300 350)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot local haplotype trajectories')

    args = parser.parse_args()
    pname = args.patient
    roi = args.roi
    VERBOSE = args.verbose
    use_plot = args.plot

    patient = load_patient(pname)
    patient.discard_nonsequenced_samples()

    # TODO: some more work if the fragment is "genomewide" or "HXB2" or so
    (fragment, start, end) = get_fragmented_roi(patient, roi, VERBOSE=VERBOSE)
    # FIXME
    sys.exit()


    refseq = patient.get_reference(fragment)
    refroi = ''.join(refseq[start: end])

    if VERBOSE >= 1:
        print patient.name, fragment, start, end

    if VERBOSE >= 2:
        print 'Get haplotype trajectories'
    (ht, indt, htseqs) = patient.get_local_haplotype_count_trajectories(\
                                        fragment, start, end, VERBOSE=VERBOSE)

    if VERBOSE >= 2:
        print 'Align haplotypes and delete rare ones'
    indht = (ht > 5).any(axis=0)
    ht = ht[:, indht]
    htseqs = htseqs[indht]
    alim = np.array(build_msa(htseqs))
    alibin = alim != alim[ht[0].argmax()]

    if VERBOSE >= 2:
        print 'Restrict to polymorphic sites'
    # FIXME: use a better criterion
    ind_poly = (alibin.sum(axis=0) > 5)
    y = np.array(alibin[:, ind_poly], float)

    # Weight with the frequencies (from all time points??) and subtract mean
    # FIXME: decide on these weights!
    # for now normalize each time point by the coverange and weight equally
    hft = (1.0 * ht.T / ht.sum(axis=1)).T
    w = (hft.sum(axis=0))**(0.3)
    yw = (w * y.T).T
    ym = yw - yw.mean(axis=0)

    if VERBOSE >= 2:
        print 'Diagonalize on covariance matrix'
    # NOTE: V[i] is the i-th eigenvector of the covariance matrix y.T x y, i.e.
    # V[i, j] is the projection of site j onto eigenvector i. For a sequence s,
    # its projection onto V[i] is sum_j (s[j] * V[i, j]), i.e. s V[i]. This is
    # also equal to U[i].
    U, s, V = linalg.svd(ym, full_matrices=False)

    if use_plot:
        # Divide the points by the weight of the haplotype, to fan them out a bit
        M3 = U[:, :3].T# / w
        tis = ht.argmax(axis=0)
        colors = cm.jet(1.0 * ht.argmax(axis=0) / len(indt))

        # 2D plot
        fig, ax = plt.subplots()
        ax.scatter(M3[0], M3[1], color=colors)
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_title(', '.join([patient.code] + map(str, roi)))
        ax.grid(True)

        ## 3D plot
        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        #ax.scatter(M3[0], M3[1], M3[2], color=colors)
        #ax.set_xlabel('PC1')
        #ax.set_ylabel('PC2')
        #ax.set_zlabel('PC3')
        #ax.set_title(', '.join([patient.code] + map(str, roi)))

        # Haplotype trajectories
        from hivwholeseq.patients.get_local_haplotypes import plot_haplotype_frequencies
        plot_haplotype_frequencies(patient.times[indt], hft,
                                   title=patient.name+', '+' '.join(map(str, roi)))

        plt.ion()
        plt.show()

