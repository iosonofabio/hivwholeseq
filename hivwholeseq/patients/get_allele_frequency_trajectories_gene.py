#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/11/13
content:    Collect the allele frequencies of all samples to the initial
            consensus and save them as a matrix into a single trajectory file.
'''
# Modules
import os
import argparse
from operator import attrgetter, itemgetter
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio import SeqIO

from hivwholeseq.miseq import alpha
from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file, \
        get_allele_counts_insertions_from_file_unfiltered, \
        filter_nus
from hivwholeseq.patients.patients import load_patient, convert_date_deltas_to_float
from hivwholeseq.patients.filenames import get_initial_reference_filename
from hivwholeseq.patients.one_site_statistics import plot_allele_frequency_trajectories_from_counts as plot_nus_from_act
from hivwholeseq.patients.one_site_statistics import plot_allele_frequency_trajectories_from_counts_3d as plot_nus_from_act_3d
from hivwholeseq.patients.one_site_statistics import get_allele_count_trajectories
from hivwholeseq.genome_info import genes, locate_gene, gene_edges
from hivwholeseq.generic_utils import mkdirs
from hivwholeseq.sequencing.primer_info import fragments_genes, \
        fragments_RNA_structures, fragments_other
from hivwholeseq.sequence_utils import correct_genbank_features_load



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Update initial consensus')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--features', required=True, nargs='+',
                        help='Feature to analyze (e.g. gag, RRE)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the allele frequency trajectories')
    parser.add_argument('--logit', action='store_true',
                        help='use logit scale (log(x/(1-x)) in the plots')
    parser.add_argument('--PCR1', type=int, default=1,
                        help='Take only PCR1 samples [0=off, 1=where both available, 2=always]')
    parser.add_argument('--threshold', type=float, default=0.9,
                        help='Minimal frequency plotted')
    
    args = parser.parse_args()
    pname = args.patient
    featnames = args.features
    VERBOSE = args.verbose
    plot = args.plot
    use_logit = args.logit
    use_PCR1 = args.PCR1
    threshold = args.threshold

    patient = load_patient(pname)
    patient.discard_nonsequenced_samples()
    refseq = SeqIO.read(patient.get_reference_filename('genomewide', 'gb'), 'gb')
    correct_genbank_features_load(refseq)
    samplenames = patient.samples.index

    fragments_features = {}
    fragments_features.update(fragments_genes)
    fragments_features.update(fragments_RNA_structures)
    fragments_features.update(fragments_other)
    fragments_features = {featname: frags
                          for (featname, frags) in fragments_features.iteritems()
                          if featname in featnames}

    for featname, fragments in fragments_features.iteritems():
        # FIX tat/rev they do not belong to a single fragment continuously
        if featname in ('tat', 'rev'):
            continue

        for fragment in fragments:
            if VERBOSE >= 1:
                print featname, fragment,

            # Collect allele counts from patient samples, and return only positive hits
            # sns contains sample names and PCR types
            (sns, act) = get_allele_count_trajectories(pname, samplenames, fragment,
                                                       use_PCR1=use_PCR1, VERBOSE=VERBOSE)
            ind = [i for i, (_, sample) in enumerate(patient.samples.iterrows())
                   if sample.name in map(itemgetter(0), sns)]
            samples = patient.samples.iloc[ind]
            times = convert_date_deltas_to_float(samples.date - patient.transmission_date,
                                                 unit='day')
            ntemplates = samples['n templates']

            # Restrict to feature
            coord_frag = None
            coord_feature = None
            for feature in refseq.features:
                if feature.id == featname:
                    coord_feature = feature.location
                if feature.id == fragment:
                    coord_frag = feature.location
                if None not in (coord_frag, coord_feature):
                    break
            else:
                continue

            frag_start = coord_frag.nofuzzy_start
            feat_start = max(0, coord_feature.nofuzzy_start - frag_start)
            feat_end = min(act.shape[2], coord_feature.nofuzzy_end - frag_start)
            act = act[:, :, feat_start: feat_end]
            if VERBOSE >= 1:
                print feat_start, feat_end

            # FIXME: use masked arrays?
            aft = (1.0 * act.swapaxes(0, 1) / act.sum(axis=1)).swapaxes(0, 1)
            aft[np.isnan(aft)] = 0
            aft[(aft < 1e-4)] = 0

            if plot is not None:
                import matplotlib.pyplot as plt

                title='Patient '+pname+', '+featname+' ('+fragment+')'
                options = []
                if featname in genes:
                    options.append('syn-nonsyn')
                if plot in ('2D', '2d', ''):

                    plot_nus_from_act(times, act,
                                      title=title,
                                      VERBOSE=VERBOSE, logit=use_logit,
                                      ntemplates=ntemplates,
                                      threshold=threshold,
                                      options=options)

                elif plot in ('3D', '3d'):
                    plot_nus_from_act_3d(times, act,
                                         title=title,
                                         VERBOSE=VERBOSE, logit=use_logit,
                                         threshold=threshold,
                                         options=options)

        if plot is not None:
            plt.tight_layout()
            plt.ion()
            plt.show()
