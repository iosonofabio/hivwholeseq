# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/11/13
content:    After the preliminary mapping to reference, plot coverage and allele
            frequencies to spot major issues.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib import cm

from hivwholeseq.miseq import read_types
from hivwholeseq.filenames import get_premapped_filename, get_reference_premap_filename, \
        get_fragment_positions_filename
from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file_unfiltered
from hivwholeseq.primer_info import primers_coordinates_HXB2_inner as pcis_HXB2
from hivwholeseq.primer_info import primers_coordinates_HXB2_outer as pcos_HXB2
from hivwholeseq.mapping_utils import get_number_reads

from hivwholeseq.samples import load_sequencing_run, SampleSeq



# Functions
def plot_coverage(counts, frags_pos=None, frags_pos_out=None,
                  title=None):
    '''Plot the coverage and the minor allele frequency'''
    cov = counts.sum(axis=1)
    cov_tot = cov.sum(axis=0)

    fig, ax = plt.subplots(1, 1, figsize=(11, 8))
    ax.plot(cov_tot.T, lw=2, c='k', label=read_types)
    ax.set_xlabel('Position [bases]')
    ax.set_ylabel('Coverage')

    # If the fragments positions are marked, plot them
    # Inner primers
    if frags_pos is not None:
        for i, frag_pos in enumerate(frags_pos.T):
            ax.plot(frag_pos, 2 * [(0.97 - 0.03 * (i % 2)) * ax.get_ylim()[1]],
                    c=cm.jet(int(255.0 * i / len(frags_pos.T))), lw=2)

    # Outer primers
    if frags_pos_out is not None:
        for i, frag_pos in enumerate(frags_pos_out.T):
            ax.plot(frag_pos, 2 * [(0.96 - 0.03 * (i % 2)) * ax.get_ylim()[1]],
                    c=cm.jet(int(255.0 * i / len(frags_pos_out.T))), lw=2)

    ax.set_xlim(-500, 9500)

    if title is not None:
        ax.set_title(title, fontsize=18)

    plt.tight_layout(rect=(0, 0, 1, 0.95)) 


def check_premap(data_folder, adaID, fragments, seq_run, samplename,
                 qual_min=30, match_len_min=10,
                 maxreads=-1, VERBOSE=0,
                 title=None):
    '''Check premap to reference: coverage, etc.'''
    refseq = SeqIO.read(get_reference_premap_filename(data_folder, adaID), 'fasta')

    fragpos_filename = get_fragment_positions_filename(data_folder, adaID)
    if os.path.isfile(fragpos_filename):
        # Load the fragment positions, considering mixed fragments (e.g. F5a+b)
        fragtmp = []
        postmp = []
        with open(fragpos_filename, 'r') as f:
            f.readline() #HEADER
            for line in f:
                fields = line[:-1].split('\t')
                fragtmp.append(fields[0])
                if 'inner' not in fields[1]:
                    postmp.append([fields[1], fields[4]])
                else:
                    start = int(fields[1].split(',')[1].split(': ')[1].rstrip('}'))
                    end = int(fields[4].split(',')[1].split(': ')[1].rstrip('}'))
                    postmp.append([start, end])

        postmp = np.array(postmp, int)

        frags_pos = np.array([postmp[fragtmp.index(fr)] for fr in fragments], int).T

    else:
        frags_pos = None
    
    frags_pos_out = None

    # Open BAM and scan reads
    input_filename = get_premapped_filename(data_folder, adaID, type='bam')

    # Count reads if requested
    if VERBOSE:
        print 'N. of reads:', get_number_reads(input_filename)

    # Get counts
    counts, inserts = get_allele_counts_insertions_from_file_unfiltered(input_filename,
                                                             len(refseq),
                                                             qual_min=qual_min,
                                                             match_len_min=match_len_min,
                                                             maxreads=maxreads,
                                                             VERBOSE=VERBOSE)

    # Plot results
    if title is None:
        title=', '.join(map(lambda x: ' '.join([x[0], str(x[1])]),
                            [['run', seq_run],
                             ['adaID', adaID],
                             ['sample', samplename],
                             ['n_reads', maxreads],
                            ]))
    plot_coverage(counts,
                  frags_pos=frags_pos,
                  frags_pos_out=frags_pos_out,
                  title=title)

    # SAVEFIG
    from hivwholeseq.adapter_info import foldername_adapter
    plt.savefig(data_folder+foldername_adapter(adaID)+'figures/coverage_premapped_'+samplename+'.png')

    return (counts, inserts)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Check consensus')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='+',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--maxreads', type=int, default=1000,
                        help='Number of reads analyzed')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--titles', nargs='*', default=None,
                        help='Give a title to the figure')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    maxreads = args.maxreads
    VERBOSE = args.verbose
    titles = args.titles

    # Specify the dataset
    dataset = load_sequencing_run(seq_run)
    data_folder = dataset.folder

    # If the script is called with no adaID, iterate over all
    samples = dataset.samples
    if adaIDs is not None:
        samples = samples.loc[samples.adapter.isin(adaIDs)]

    if VERBOSE >= 3:
        print 'adaIDs', samples.adapter

    for i, (samplename, sample) in enumerate(samples.iterrows()):
        sample = SampleSeq(sample)
        adaID = sample.adapter
        fragments = sample.regions_complete

        if VERBOSE:
            print seq_run, adaID
        if titles is not None:
            title = titles[i]
        else:
            title = None

        (counts, inserts) = check_premap(data_folder, adaID,
                                         fragments, seq_run, samplename,
                                         maxreads=maxreads,
                                         VERBOSE=VERBOSE,
                                         title=title)

    plt.ion()
    plt.show()

