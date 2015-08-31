# vim: fdm=marker
'''
author:     Fabio Zanini
date:       19/05/14
content:    Calculate and plot the propagator of allele frequencies.
'''
# Modules
import sys
import argparse
from operator import attrgetter
import numpy as np
from Bio import SeqIO
from Bio.Seq import MutableSeq, translate
from Bio.Alphabet.IUPAC import unambiguous_dna
import matplotlib.pyplot as plt
from seqanpy import align_overlap

from hivwholeseq.utils.miseq import alpha
from hivwholeseq.reference import load_custom_reference
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.propagator_allele_frequency import Propagator



# Functions
def get_gene_positions_in_fragment(gene, fragment, gwseq, fragseq, VERBOSE=0):
    '''Get the coordinates of a gene within a fragment'''
    # Find coordinates of gene in reference
    feagene = gwseq.features[map(attrgetter('id'), gwseq.features).index(gene)]
    gene_start = feagene.location.nofuzzy_start
    gene_end = feagene.location.nofuzzy_end

    # Sanity check on coordinates
    feafrag = refseq.features[map(attrgetter('id'), gwseq.features).index(fragment)]
    fragrefgw = feafrag.extract(gwseq)
    if len(fragseq) != len(fragrefgw):
        raise ValueError('Problem with coordinates between fragment and genomewide.')

    # Find coordinates of gene in fragment
    frag_start = feafrag.location.nofuzzy_start
    frag_end = feafrag.location.nofuzzy_end

    # complete gene
    if (frag_start <= gene_start) and (frag_end >= gene_end):
        if VERBOSE >= 2:
            print 'Complete gene found'
        positions = np.arange(gene_start, gene_end) - frag_start

    # start of gene
    elif (frag_start <= gene_start):
        if VERBOSE >= 2:
            print 'WARNING: only gene start found'
        positions = np.arange(gene_start, frag_end) - frag_start
        if len(positions) % 3:
            positions = positions[:-(len(positions) % 3)]

    # end of gene
    elif (frag_end >= gene_end):
        if VERBOSE >= 2:
            print 'WARNING: only gene end found'
        positions = np.arange(frag_start, gene_end) - frag_start
        if len(positions) % 3:
            positions = positions[len(positions) % 3:]

    # middle of gene: guess reading frame
    else:
        if VERBOSE >= 2:
            print 'WARNING: only gene middle found'
        prot = feagene.extract(gwseq).seq.translate()
        ali_score = []
        for rf_start in xrange(3):
            tmpseq = fragseq[rf_start:].seq
            if len(tmpseq) % 3:
                tmpseq = tmpseq[:-(len(tmpseq) % 3)]
            tmpprot = tmpseq.translate()
            (score, ali1, ali2) = align_overlap(tmpprot, prot)
            ali_score.append(score)
        rf_start = np.argmax(ali_score)
        positions = np.arange(frag_start + rf_start, frag_end) - frag_start
        if len(positions) % 3:
            positions = positions[:-(len(positions) % 3)]

    return positions



# Script
if __name__ == '__main__': 

    parser = argparse.ArgumentParser(description='Propagator for allele frequencies',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save the propagator to file')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the propagator')
    parser.add_argument('--deltat', type=int, nargs=2, default=[100, 300],
                        help='Time in days between final and initial (range)')
    parser.add_argument('--logit', action='store_true',
                        help='use logit scale (log(x/(1-x)) in the plots')
    parser.add_argument('--min-depth', type=int, default=100, dest='min_depth',
                        help='Minimal depth to consider the site')
    parser.add_argument('--genes', default=('gag', 'pol'), nargs='+',
                        help='Gene to analyze (e.g. pol)')

    args = parser.parse_args()
    genes = map(lambda x: x.lower(), args.genes)
    pnames = args.patients
    VERBOSE = args.verbose
    save_to_file = args.save
    plot = args.plot
    dt = args.deltat
    use_logit = args.logit
    depth_min = args.min_depth

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]
        # FIXME: the initial ref of p7 is mislabelled and a mess
    else:
        patients = patients.loc[patients.code != 'p7']

    # Prepare output structures
    n_binsx = 5
    binsy = [0.,
             0.002,
             0.01,
             0.025,
             0.12,
             0.33,
             0.67,
             0.88,
             0.975,
             0.99,
             0.998,
             1.]
    props = {(gene, synkey): Propagator(n_binsx, binsy=binsy, use_logit=use_logit)
             for gene in genes for synkey in ('syn', 'nonsyn')}

    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        samplenames = patient.samples.index

        refseq = patient.get_reference('genomewide', format='gb')

        for gene in genes:

            if VERBOSE >= 1:
                print pname, gene,

            # Get the right fragment(s)
            # FIXME: do better than this ;-)
            frags = {'pol': ['F2', 'F3'], 'gag': ['F1'], 'env': ['F5', 'F6']}
            fragments = frags[gene]

            if VERBOSE >= 1:
                print fragments

            for fragment in fragments:
                fragseq = patient.get_reference(fragment)
                # FIXME: there are sometimes coordinate shifts in the genomewide ref.
                try:
                    positions = get_gene_positions_in_fragment(gene, fragment, refseq, fragseq,
                                                               VERBOSE=VERBOSE)
                except ValueError:
                    continue

                # Get allele freqs
                aft, ind = patient.get_allele_frequency_trajectories(fragment,
                                                                     cov_min=depth_min)
                aft = aft[:, :, positions]

                n_templates = np.array(patient.n_templates[ind])
                indd = n_templates >= depth_min
                aft = aft[indd]
                ind = ind[indd]
                n_templates = n_templates[indd]

                ts = patient.times[ind]

                # Make initial consensus to declare syn/nonsyn
                # NOTE: the initial consensus is NOT the mapping reference, so there
                # might be gaps in the sequence. That's ok.
                genem = alpha[aft[0].argmax(axis=0)]

                # Collect counts
                for i in xrange(aft.shape[0] - 1):
                    for j in xrange(i + 1, aft.shape[0]):
                        dij = j - i
                        if ((ts[j] - ts[i]) > dt[1]) or ((ts[j] - ts[i]) < dt[0]):
                            continue

                        afaf = {'syn': [], 'nonsyn': []}

                        for pos in xrange(aft.shape[2]):
                            pos_cod = pos // 3
                            pos_wn = pos % 3

                            # Because the initial consensus is NOT the mapping reference,
                            # genem might have gaps... that's ok, they come in groups of three.
                            # FIXME: this is slightly more messed up than that, for now skip the
                            # whole codon
                            if '-' in genem[pos_cod * 3: (pos_cod + 1) * 3]:
                                continue

                            # TODO: Check this algorithm!!
                            for ia, a in enumerate(alpha[:4]):

                                # exclude ancestral alleles
                                if a == genem[pos]:
                                    continue

                                cod_anc = genem[pos_cod * 3: (pos_cod + 1) * 3]
                                cod_new = cod_anc.copy()
                                cod_new[pos_wn] = a

                                if translate(''.join(cod_anc)) == translate(''.join(cod_new)):
                                # FIXME: random control
                                #if np.random.rand() > 0.5:
                                    afaf['syn'].append((aft[i, ia, pos], aft[j, ia, pos]))
                                else:
                                    afaf['nonsyn'].append((aft[i, ia, pos], aft[j, ia, pos]))

                        for key in ('syn', 'nonsyn'):
                            pp = props[(gene, key)]
                            pp.histogram += np.histogram2d(*np.transpose(afaf[key]),
                                                           bins=[pp.binsx, pp.binsy])[0]


    if plot:
        fig, ax = plt.subplots(figsize=(9, 6))
        ax.set_title(gene+', $\Delta t \in ['+str(dt[0])+', '+str(dt[1])+']$ days')
        props[(gene, 'syn')].plot(heatmap=False, figaxs=(fig, ax), ls='--', marker='s')
        props[(gene, 'nonsyn')].plot(heatmap=False, figaxs=(fig, ax), start_nu=False)

        plt.ion()
        plt.show()

