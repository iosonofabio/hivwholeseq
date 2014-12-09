# vim: fdm=marker
'''
author:     Fabio Zanini
date:       09/12/14
content:    Compare frequency of minor variants in PCR1/2 for the test cell sample.
'''
# Modules
import os

from hivwholeseq.patients.samples import load_sample_sequenced
from hivwholeseq.patients.get_local_haplotypes import cluster_haplotypes



# Script
if __name__ == '__main__':

    fragment = 'F4'
    start = 420
    end = 540
    covmin = 5

    sample1 = load_sample_sequenced('20097-C1')
    sample2 = load_sample_sequenced('20097-C1-2')


    h1 = cluster_haplotypes(sample1.get_local_haplotypes(fragment, start, end, PCR=1),
                            min_abundance=covmin)
    h2 = cluster_haplotypes(sample2.get_local_haplotypes(fragment, start, end, PCR=2),
                            min_abundance=covmin)

    seqs_set = set(h1.keys()) | set(h2.keys())
    tot1 = sum(h1.values())
    tot2 = sum(h2.values())

    import hivwholeseq.plot_utils
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    for seq in seqs_set:
        ax.scatter(1.0 * h1[seq] / tot1, 1.0 * h2[seq] / tot2, s=50,
                   c='b', edgecolor='none')

    xmin = 1e-3
    ax.plot([xmin, 1 - xmin], [xmin, 1 - xmin], c='k', lw=2)

    ax.set_xlabel('Haplotype freq PCR1')
    ax.set_ylabel('Haplotype freq PCR2')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(xmin, 0.5)
    ax.set_ylim(xmin, 0.5)
    ax.grid(True)
    ax.set_title('Cell sample test, PCR1 vs PCR2')

    plt.tight_layout()

    plt.ion()
    plt.show()
