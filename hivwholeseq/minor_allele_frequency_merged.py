# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/12/13
content:    Plot the minor allele frequencies, genome wide.
'''
# Modules
import argparse
import numpy as np

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.samples import samples
from hivwholeseq.filenames import get_merged_allele_frequencies_filename



# Functions
def plot_minor_allele_frequency_filtered(data_folder, adaID, fragments, VERBOSE=0,
                                savefig=False):
    '''Plot minor allele frequency along the genome''' 
    nus = np.load(get_merged_allele_frequencies_filename(data_folder, adaID, fragments))

    nu_min = np.ma.masked_all(nus.shape[-1])
    for pos, nutmp in enumerate(nus.T):
        try:
            if not np.ma.is_masked(nutmp):
                nu_min[pos] = np.sort(nutmp)[-2]
        except ValueError:
            print pos, np.ma.is_masked(nutmp)
            import ipdb; ipdb.set_trace()

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(15, 8))
    ax.plot(nu_min, lw=1.5, c='k')
    ax.scatter(np.arange(len(nu_min)), nu_min, s=30, c='k')
    ax.set_yscale('log')
    ax.set_xlabel('Position')
    ax.set_ylabel(r'$\nu$', fontsize=20)
    ax.set_title('adaID '+adaID+', '+'-'.join(fragments))
    ax.set_xlim(-100, len(nu_min) + 100)

    plt.tight_layout()

    if savefig:
        from hivwholeseq.filenames import \
                get_minor_allele_frequency_merged_figure_filename as gff
        outputfile = gff(data_folder, adaID, fragments)
        fig.savefig(outputfile)
        plt.close(fig)
    else:
        plt.ion()
        plt.show()



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Study minor allele frequency')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--no-savefig', action='store_false', dest='savefig',
                        help='Show figure instead of saving it')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose
    savefig = args.savefig

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = dataset['adapters']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # Iterate over all requested samples
    for adaID in adaIDs:

        # If the script is called with no fragment, iterate over all
        samplename = dataset['samples'][dataset['adapters'].index(adaID)]
        if not fragments:
            fragments_sample = [fr[:2] for fr in samples[samplename]['fragments']]
        else:
            from re import findall
            fragments_all = samples[samplename]['fragments']
            fragments_sample = []
            for fragment in fragments.split('-'):
                frs = filter(lambda x: fragment in x, fragments_all)
                if len(frs):
                    fragments_sample.append(frs[0][:2])

        plot_minor_allele_frequency_filtered(data_folder, adaID, fragments_sample,
                                             VERBOSE=VERBOSE, savefig=savefig)
