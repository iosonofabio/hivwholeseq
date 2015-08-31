# vim: fdm=marker
'''
author:     Fabio Zanini
date:       31/03/14
content:    Plot the filtered minor allele frequency for a clone and a patient
            sample.
'''
# Modules
import argparse
from operator import itemgetter
from itertools import izip
import numpy as np
import matplotlib.pyplot as plt

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.utils.miseq import alpha, read_types
from hivwholeseq.sequencing.filenames import get_allele_frequencies_filename




# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Study minor allele frequency')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    fragments = args.fragments
    VERBOSE = args.verbose

    # Reference sample
    d_ref = {'seqrun': 'Tue28', 'adaID': 'TS2'}
    d_pat = {'seqrun': 'Tue48', 'adaID': 'N2-S3'}
    d_both = [d_ref, d_pat]

    for fragment in fragments:
        fig, axs = plt.subplots(1, 2, figsize=(14, 8))

        for ip, (d, ax) in enumerate(izip(d_both, axs)):
            dataset = MiSeq_runs[d['seqrun']]
            data_folder = dataset['folder']
            table = dataset['sample_table']
            adaID = d['adaID']
            name = table.set_index('adaID', drop=False).loc[adaID]['name']

            nu_filtered = np.load(get_allele_frequencies_filename(data_folder, adaID, fragment))
            nut = np.zeros(nu_filtered.shape[-1])
            for pos, nupos in enumerate(nu_filtered.T):
                nut[pos] = np.sort(nupos)[-2]

            ax.set_yscale('log')
            ax.set_title(name, fontsize=16)
            ax.set_xlabel('Position in HIV region [bp]', fontsize=16)
            if ip == 0:
                ax.set_ylabel('Frequency of minor allele', fontsize=16)
            
            # Plot filtered
            ax.plot(nut, lw=1.5, c='b',
                    alpha=0.7, label='Filtered')
            ax.scatter(np.arange(len(nut)),
                       nut, lw=1.5, c='b',
                       alpha=0.7)

            ax.set_xlim(-100, len(nut) + 100)
            ax.set_ylim(1e-5, 1)
        
        plt.tight_layout(rect=(0, 0, 1, 0.99))
 
