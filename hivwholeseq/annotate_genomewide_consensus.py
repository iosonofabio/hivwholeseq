# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/12/13
content:    Add gene and region annotations to the genome wide consensus.
'''
# Modules
import argparse
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.samples import samples
from hivwholeseq.filenames import get_merged_consensus_filename as gmcf


# Functions
def extract_feature(refseq, featurename):
    '''Extract a feature from a SeqRecord'''
    if not len(refseq.features):
        raise ValueError('The sequence has no features')

    from operator import attrgetter
    featurenames = map(attrgetter('id'), refseq.features)
    try:
        record = refseq.features[featurenames.index(featurename)].extract(refseq)
        record.id = record.id+'|'+featurename
        record.name = record.name+'|'+featurename
    except ValueError:
        raise ValueError('Feature not in the list of features of the sequence')

    return record


def annotate_sequence(seqrecord, features=['gene', 'RNA structure', 'other']):
    '''Annotate a consensus with the genes and stuff (in place)'''
    from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
    from hivwholeseq.genome_info import gene_edges, RNA_structure_edges, \
            other_edges, find_region_edges, find_region_edges_multiple
    from hivwholeseq.primer_info import primers_PCR as primers_PCR_edges
    edge_dict = {'gene': gene_edges,
                 'RNA structure': RNA_structure_edges,
                 'PCR primers': primers_PCR_edges,
                 'other': other_edges}

    smat = np.array(seqrecord)

    for feature_type in features:
        edges_all = edge_dict[feature_type]
        for name, edges in edges_all.iteritems():
            # Skip a feature if it's present already
            if name in map(lambda x: x.id, seqrecord.features):
                continue

            # Behave differently for unsplit regions and split ones
            if len(edges) == 2:
                # LTR problems with F6
                if 'F6' in name:
                    pos_edge = find_region_edges(smat[::-1], [edges[1][::-1], edges[0][::-1]])
                    pos_edge = [len(smat) - 1 - pos_edge[1], len(smat) - 1 - pos_edge[0]]
                else:
                    pos_edge = find_region_edges(smat, edges)
                location = FeatureLocation(*pos_edge)
            else:
                pos_edges = find_region_edges_multiple(smat, edges)
                locations = [FeatureLocation(*pos_edge) for pos_edge in pos_edges]
                location = CompoundLocation(locations)
            feature = SeqFeature(location, type=feature_type, id=name, strand=1)
            seqrecord.features.append(feature)


def annotate_plot(ax, consensus, features=['gene', 'RNA structure', 'other']):
    '''Annotate a plot with genome features'''
    from Bio.SeqFeature import CompoundLocation
    from matplotlib.patches import Rectangle
    from matplotlib import cm

    annotate_sequence(consensus)
    feats = filter(lambda x: x.type in features, consensus.features)

    (ymin, ymax) = ax.get_ylim()
    if ax.get_yscale() == 'linear':
        yspan = ymax - ymin
        yannos = [ymin - (0.05 + 0.03 * i) * yspan for i in xrange(len(feats))]
        ymin_new = ymin - 0.4 * yspan
    else:
        yspan = ymax / ymin
        yannos = [np.exp(np.log(ymin) - (0.05 + 0.03 * i) * np.log(yspan))
                  for i in xrange(len(feats))]
        ymin_new = ymin / np.exp(0.4 * np.log(yspan))

    for i, feat in enumerate(feats):

        if not isinstance(feat.location, CompoundLocation):
            xstart, xend = feat.location.start, feat.location.end
            xspan = xend - xstart
            ax.add_patch(Rectangle((xstart, ymin_new), xspan, ymax - ymin_new,
                                   fc=cm.jet(int(255.0 * i / len(feats))),
                                   alpha=0.5,
                                   ec='none'))
            ax.plot([xstart, xend], 2 * [yannos[i]], lw=2.5,
                     c='k', alpha=0.5)
            ax.text(xstart + 0.1 * xspan, yannos[i], feat.id)

        else:
            for j, location in enumerate(feat.location.parts):
                xstart, xend = location.start, location.end
                xspan = xend - xstart
                ax.add_patch(Rectangle((xstart, ymin_new), xspan, ymax - ymin_new,
                                       fc=cm.jet(int(255.0 * i / len(feats))),
                                       alpha=0.5,
                                       ec='none'))
                ax.plot([xstart, xend], 2 * [yannos[i]], lw=2.5,
                         c='k', alpha=0.5)
                ax.text(xstart + 0.1 * xspan, yannos[i], feat.id+str(j+1))

    ax.set_ylim(ymin_new, ymax)

            

# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Merge consensi and allele frequencies')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    
    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    VERBOSE = args.verbose

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = dataset['adapters']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    for adaID in adaIDs:

        samplename = dataset['samples'][dataset['adapters'].index(adaID)]
        fragments = [fr[:2] for fr in samples[samplename]['fragments']]

        # Annotate sequence
        consall = SeqIO.read(gmcf(data_folder, adaID, fragments), 'fasta')

        # Plot minor allele frequencies and annotate plot
        from hivwholeseq.minor_allele_frequency_merged import \
                plot_minor_allele_frequency_filtered
        plot_minor_allele_frequency_filtered(data_folder, adaID, fragments,
                                             VERBOSE=VERBOSE, savefig=False)

        ax = plt.gca()

        annotate_plot(ax,  consall)

        plt.ion()
        plt.show()



