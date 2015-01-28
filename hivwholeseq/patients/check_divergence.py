# vim: fdm=marker
'''
author:     Fabio Zanini
date:       26/03/14
content:    Check divergence of sequences during time in various ways.
'''
# Modules
import os
import argparse
import numpy as np
import pysam
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.patients import load_patient
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_mapped_to_initial_filename, get_allele_frequency_trajectories_filename, \
        get_allele_count_trajectories_filename
from hivwholeseq.utils.mapping import extract_mapped_reads_subsample_open, pair_generator


# Functions
def get_distance_reads_sequence(seq, reads, VERBOSE=0,
                                score_match=3,
                                score_mismatch=-3):
    '''Get the distance in alignment score between read pairs and a sequence'''
    from seqanpy import align_overlap

    seqs = ''.join(seq)
    deltas = []
    for irp, read_pair in enumerate(reads):
        d = 0
        for read in read_pair:
            (score, alis, alir) = align_overlap(seqs, read.seq,
                                                score_match=score_match,
                                                score_mismatch=score_mismatch)
            start = len(alir) - len(alir.lstrip('-'))
            end = len(alir.rstrip('-'))
            scoremax = score_match * (end - start)
            delta = scoremax - score
            d += delta

        deltas.append(d)

    return deltas



# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Check patient samples')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--allele-frequencies', action='store_true', dest='afs',
                        help='Use allele frequencies')
    parser.add_argument('--reads', type=int, default=0,
                        help='Use reads')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose
    use_af = args.afs
    maxreads = args.reads

    patient = load_patient(pname)
    times = patient.times

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    if VERBOSE:
        print pname

    if use_af:
        divs = []
        for fragment in fragments:

            aft_filename = get_allele_frequency_trajectories_filename(pname, fragment)
            aft = np.ma.array(np.load(aft_filename))
            aft[aft < 1e-2] = np.ma.masked

            cons = alpha[aft[0].argmax(axis=0)]
            anc = frozenset(zip(cons, np.arange(len(cons))))

            from itertools import izip
            n_muts = []
            for af in aft:
                muts = frozenset([(alpha[ai], pos)
                                  for (ai, pos) in izip(*(af > 0.5).nonzero())
                                  if alpha[ai] != '-'])
                muts -= anc
                n_muts.append(len(muts))

            div = (0.5 * n_muts[-1] / times[-1] / len(cons))
            div += (0.5 * n_muts[-2] / times[-2] / len(cons))
            divs.append(div)
            if VERBOSE:
                print fragment, '{:1.2e}'.format(div * 365.25), 'substitutions per site per year'

            # Divergence at low frequency
            nutots = []
            for af in aft:
                nutot = sum(1 for (ai, pos) in izip(*((af > 0.01) & (af < 0.1)).nonzero())
                            if (alpha[ai] != '-') and (cons[pos] != alpha[ai]))
                nutots.append(nutot)
            print fragment, 'nutots', nutots 

    if maxreads != 0:
        for fragment in fragments:
            refseq = patient.get_reference(fragment)
            dists = {}
            for sample in patient.itersamples():
                if VERBOSE >= 1:
                    print sample.name,

                if sample[fragment] not in ['ok', 'low']:
                    if VERBOSE >= 1:
                        print 'not "ok". skipping'
                    continue

                if VERBOSE >= 1:
                    print 'ok'
                               
                bamfilename = sample.get_mapped_filtered_filename(fragment, decontaminated=True)
                with pysam.Samfile(bamfilename, 'rb') as bamfile:
                    if maxreads == -1:
                        reads = pair_generator(bamfile)
                    else:
                        reads = extract_mapped_reads_subsample_open(bamfile, maxreads,
                                                                    VERBOSE=VERBOSE)

                    dists[sample.name] = get_distance_reads_sequence(refseq, reads,
                                                                     VERBOSE=VERBOSE,
                                                                     score_match=3,
                                                                     score_mismatch=-3)


            hs = {}
            binmax = max(map(max, dists.itervalues()))
            bins = np.arange(0, binmax, 6)
            bincenters = 0.5 * (bins[1:] + bins[:-1])
            for samplename, dist in dists.iteritems():
                hs[samplename] = np.histogram(dist, bins=bins, density=True)[0]

            fig, ax = plt.subplots(figsize=(10, 5))
            ax.set_xlabel('Divergence [delta alignment score]')
            ax.set_ylabel('Prob density of read pairs')
            for i, sample in enumerate(patient.itersamples()):
                if sample.name in hs:
                    h = hs[sample.name]
                    label = str(sample['days since infection'])+', '+sample.name
                    ax.plot(bincenters, h, lw=2, label=label,
                            color=cm.jet(1.0 * i / len(patient.samples)),
                            alpha=0.8)

            ax.grid(True)
            ax.legend(loc=1, title='Days since infection:', fontsize=14)
            ax.set_title(patient.name+', '+fragment)
            plt.tight_layout()

    plt.ion()
    plt.show()

