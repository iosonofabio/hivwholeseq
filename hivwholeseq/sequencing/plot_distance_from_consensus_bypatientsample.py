# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/10/14
content:    In some samples, the cumulative distribution of distance from consensus
            has a shoulder around 15 mutations/read pair. What are those reads?
'''
# Modules
import argparse
from operator import itemgetter
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.sequencing.check_distance_mapped_consensus import get_distance_histogram, plot_distance_histogram
from hivwholeseq.patients.samples import load_samples_sequenced, SamplePat
from hivwholeseq.sequencing.samples import SampleSeq
from hivwholeseq.generic_utils import mkdirs



# Globals



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Plot dist from consensus distribution',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    pats_or_samples = parser.add_mutually_exclusive_group(required=True)
    pats_or_samples.add_argument('--patients', nargs='+',
                                 help='Patient to analyze')
    pats_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to map (e.g. VL98-1253 VK03-4298)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save the allele counts to file')

    args = parser.parse_args()
    pnames = args.patients
    samplenames = args.samples
    VERBOSE = args.verbose
    use_save = args.save

    fragments = ['F'+str(i+1) for i in xrange(6)]

    samples = load_samples_sequenced()
    if pnames is not None:
        samples = samples.loc[samples.patient.isin(pnames)]
    elif samplenames is not None:
        samples = samples.loc[samples.index.isin(samplenames)]

    for samplename, sample in samples.iterrows():
        sample = SamplePat(sample)
        if VERBOSE >= 1:
            print samplename

        dist_hists = []
        samples_seq = sample.get_sequenced_samples()
        samples_seq = samples_seq.loc[samples_seq.PCR == 1]
        for samplename_seq, sample_seq in samples_seq.iterrows():
            sample_seq = SampleSeq(sample_seq)
            data_folder = sample_seq.seqrun_folder
            adaID = sample_seq.adapter

            for fragment in fragments:
                try:
                    dist_hist = get_distance_histogram(data_folder, adaID, fragment,
                                                       VERBOSE=VERBOSE)
                except IOError:
                    continue
                dist_hists.append((samplename_seq, fragment, dist_hist))

        dist_hists.sort(key=itemgetter(1))

        fig, ax = plt.subplots()
        for i, (samplename_seq, fragment, h) in enumerate(dist_hists):
            plot_distance_histogram(h, ax=ax,
                                    color=cm.jet(1.0 * i / len(dist_hists)),
                                    label=', '.join([samplename_seq, fragment]))
        ax.set_title(samplename)
        ax.legend(loc=1, fontsize=10)

        if use_save:
            foldername = sample.get_foldername()+'figures/'
            mkdirs(foldername)
            fn = foldername+'distance_to_consensus_seqsamples.png'
            fig.savefig(fn)
            plt.close(fig)

    if not use_save:
        plt.ion()
        plt.show()
