# vim: fdm=marker
'''
author:     Fabio Zanini
date:       24/09/14
content:    Estimate cross contamination between the samples from the decontamination
            script output.
'''
# Modules
import sys
import os
import argparse
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.patients.patients import load_samples_sequenced, SamplePat
from hivwholeseq.patients.filenames import get_decontaminate_summary_filename, \
        get_crosscontamination_figure_filename
from hivwholeseq.patients.decontaminate_reads import get_number_reads_summary, refnames


# Script
if __name__ == '__main__':

    
    # Parse input args
    parser = argparse.ArgumentParser(description='Get contamination matrix',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--fragments', nargs='+',
                        default=['F'+str(i+1) for i in xrange(6)],
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', action='store_true',
                        help='Plot the matrix')
    parser.add_argument('--saveplot', action='store_true',
                        help='Save the plot')
    parser.add_argument('--compress', action='store_true',
                        help='Delete samples without a summary file')

    args = parser.parse_args()
    fragments = args.fragments
    VERBOSE = args.verbose
    plot = args.plot
    use_saveplot = args.saveplot
    use_compress = args.compress

    samples = load_samples_sequenced()
    samplenames = refnames + samples.index.tolist()

    for fragment in fragments:

        mat = np.ma.masked_all((len(samplenames), len(samplenames)))

        for (samplename, sample) in samples.iterrows():
            isp = samplenames.index(samplename)
            sample = SamplePat(sample)
            fn = get_decontaminate_summary_filename(sample.patient, samplename, fragment,
                                                    PCR=1)
            if not os.path.isfile(fn):
                continue

            n_good, n_cont, n_cont_dict = get_number_reads_summary(fn, details=True)
            
            mat[isp] = 0
            mat[isp, isp] = n_good
            for cont_name, cont_n_reads in n_cont_dict.iteritems():
                mat[isp, samplenames.index(cont_name)] = cont_n_reads

        samplenames_from = list(samplenames)
        if use_compress:
            indi_to = (-mat.mask.all(axis=1)).nonzero()[0]
            samplenames_to = [samplenames[i] for i in indi_to]
            mat = mat[indi_to]
        else:
            samplenames_to = list(samplenames)

        if plot:
            fig, ax = plt.subplots(figsize=(21, 14))
            ax.set_title('Contamination matrix, '+fragment)
            ax.set_xlabel('From:')
            ax.set_ylabel('To:')
            ax.set_xticks(np.arange(len(samplenames_from)))
            ax.set_xticklabels(samplenames_from, rotation=90, fontsize=6)
            ax.set_yticks(np.arange(len(samplenames_to)))
            ax.set_yticklabels(samplenames_to, fontsize=6)

            z = (mat.T / (mat.sum(axis=1))).T
            z_bg = 0.1 * z[z > 0].min()
            z = np.log10(z + z_bg)
            h = ax.imshow(z, cmap=cm.jet, interpolation='nearest')

            fig.colorbar(h, shrink=0.7)
            plt.tight_layout()

            if use_saveplot:
                fn_fig = get_crosscontamination_figure_filename(fragment)
                fig.savefig(fn_fig)

    if plot:
        plt.ion()
        plt.show()
