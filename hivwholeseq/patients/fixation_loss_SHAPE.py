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

from hivwholeseq.miseq import alpha
from hivwholeseq.reference import load_custom_reference
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.propagator_allele_frequency import Propagator
from hivwholeseq.shape import add_SHAPE_to_seqrecord



# Functions
def get_gene_overlaps(gene, gwseq, VERBOSE=0):
    '''Get overlap between a gene and any other one'''
    feagene = gwseq.features[map(attrgetter('id'), gwseq.features).index(gene)]
    gene_start = feagene.location.nofuzzy_start
    gene_end = feagene.location.nofuzzy_end

    from hivwholeseq.genome_info import genes
    overlap = np.zeros(gene_end - gene_start, bool)

    for feagene2 in gwseq.features:
        gene2 = feagene2.id
        if (gene2 == gene) or (gene2 not in genes):
            continue

        exon_edges = [(p.nofuzzy_start, p.nofuzzy_end) for p in feagene2.location.parts]
        for iex, (exon2_start, exon2_end) in enumerate(exon_edges):
            if (exon2_start >= gene_end) or (exon2_end <= gene_start):
                continue

            if VERBOSE >= 2:
                if len(exon_edges) == 1:
                    print gene, 'overlaps with', gene2
                else:
                    print gene, 'overlaps with', gene2, '(exon '+str(iex + 1)+')'

            ov_start = max(0, exon2_start - gene_start)
            ov_end = min(gene_end - gene_start, exon2_end - gene_start)
            overlap[ov_start: ov_end] = True

    return overlap


def get_gene_positions_in_fragment(gene, fragment, gwseq, fragseq, VERBOSE=0,
                                   overlap_info=False):
    '''Get the coordinates of a gene within a fragment'''
    # Find coordinates of gene in reference
    feagene = gwseq.features[map(attrgetter('id'), gwseq.features).index(gene)]
    gene_start = feagene.location.nofuzzy_start
    gene_end = feagene.location.nofuzzy_end

    # Sanity check on coordinates
    feafrag = refseq.features[map(attrgetter('id'), gwseq.features).index(fragment)]
    fragrefgw = feafrag.extract(gwseq)
    if len(fragseq) != len(fragrefgw):
        raise IndexError('Problem with coordinates between fragment and genomewide.')

    if overlap_info:
        overlap = get_gene_overlaps(gene, gwseq, VERBOSE=VERBOSE)

    # Find coordinates of gene in fragment
    frag_start = feafrag.location.nofuzzy_start
    frag_end = feafrag.location.nofuzzy_end

    positions = np.arange(gene_start, gene_end) - frag_start

    # complete gene
    if (frag_start <= gene_start) and (frag_end >= gene_end):
        if VERBOSE >= 2:
            print 'Complete gene found'
        ind = np.arange(len(positions))

    # start of gene
    elif (frag_start <= gene_start):
        if VERBOSE >= 2:
            print 'WARNING: only gene start found'

        ind = np.arange(frag_end - gene_start)
        if len(ind) % 3:
            ind = ind[:-(len(ind) % 3)]

    # end of gene
    elif (frag_end >= gene_end):
        if VERBOSE >= 2:
            print 'WARNING: only gene end found'

        ind = np.arange(frag_start - gene_start, gene_end - gene_start)
        if len(ind) % 3:
            ind = ind[len(ind) % 3:]

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

        ind = np.arange(frag_start + rf_start - gene_start, frag_end - gene_start)
        if len(ind) % 3:
            ind = ind[:-(len(ind) % 3)]

    positions = positions[ind]
    if overlap_info:
        overlap = overlap[ind]
        return (positions, overlap)
    else:
        return positions



# Script
if __name__ == '__main__': 

    parser = argparse.ArgumentParser(description='Fixation probability for polymorphisms',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save the results to file')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the propagator')
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
    use_logit = args.logit
    depth_min = args.min_depth

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]
        # FIXME: the initial ref of 15107 is mislabelled and a mess
    else:
        patients = patients.loc[patients.index != '15107']

    # Min-max frequencies to calculate correlations on
    xmin = 0.25
    xmax = 0.75
    xloss = 0.01
    xfix = 0.99

    # Prepare output structures
    props = {(gene, synkey, fate): []
             for gene in genes for synkey in ('syn', 'nonsyn')
             for fate in ('fix', 'loss')}

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
                    (positions, overlap) = get_gene_positions_in_fragment(gene,
                                            fragment, refseq, fragseq,
                                            VERBOSE=VERBOSE,
                                            overlap_info=True)
                except IndexError:
                    continue

                add_SHAPE_to_seqrecord(fragseq, VERBOSE=VERBOSE)
                shapefrag = fragseq.letter_annotations['SHAPE']
                if VERBOSE >= 2:
                    print 'SHAPE added'

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
                # FIXME: one should add the contraint that the mutation must BREAK
                # the pair for SHAPE to be reliable
                for pos in xrange(aft.shape[2]):
                    pos_cod = pos // 3
                    pos_wn = pos % 3

                    # Because the initial consensus is NOT the mapping reference,
                    # genem might have gaps... that's ok, they come in groups of three.
                    # FIXME: this is slightly more messed up than that, for now skip the
                    # whole codon
                    if '-' in genem[pos_cod * 3: (pos_cod + 1) * 3]:
                        continue

                    if shapefrag.mask[pos]:
                        continue

                    for ia, a in enumerate(alpha[:4]):

                        # exclude ancestral alleles
                        if a == genem[pos]:
                            continue

                        cod_anc = genem[pos_cod * 3: (pos_cod + 1) * 3]
                        cod_new = cod_anc.copy()
                        cod_new[pos_wn] = a

                        is_syn = translate(''.join(cod_anc)) == translate(''.join(cod_new))
                        is_syn &= (-overlap[pos])
                        if is_syn:
                            synkey = 'syn'
                        else:
                            synkey = 'nonsyn'

                        afp = aft[:, ia, pos]
                        ibt = ((afp[:-1] >= xmin) & (afp[:-1] <= xmax)).nonzero()[0]
                        if not len(ibt):
                            continue
                        
                        ibt = ibt[0]
                        tfix = (afp[ibt:] > xfix).nonzero()[0] + ibt
                        tlos = (afp[ibt:] < xloss).nonzero()[0] + ibt
                        if (not len(tlos)) and (not len(tfix)):
                            continue

                        elif not len(tfix):
                            fate = 'loss'

                        elif not len(tlos):
                            fate = 'fix'

                        else:
                            if tfix[0] < tlos[0]:
                                fate = 'fix'
                            else:
                                fate = 'loss'

                        props[(gene, synkey, fate)].append(shapefrag[pos])


    if plot:
        fig, ax = plt.subplots(figsize=(8, 6))
        for (gene, synkey, fate), shapes in props.iteritems():
            shapes = np.sort(shapes)
            ax.plot(shapes, np.linspace(0, 1, len(shapes)), lw=2,
                    label=', '.join([gene, synkey, fate]))

        ax.set_xlabel('SHAPE score')
        ax.set_ylabel('Cumulative distribution')
        ax.set_title('Correlation between fixation and SHAPE score')
        ax.grid(True)
        ax.set_ylim(-0.05, 1.05)
        ax.legend(loc=4, fontsize=14)

        plt.ion()
        plt.show()

