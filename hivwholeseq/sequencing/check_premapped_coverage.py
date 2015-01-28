#!/usr/bin/env python
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
from hivwholeseq.sequencing.filenames import get_premapped_filename, get_reference_premap_filename, \
        get_fragment_positions_filename
from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file_unfiltered
from hivwholeseq.sequencing.primer_info import primers_coordinates_HXB2_inner as pcis_HXB2
from hivwholeseq.sequencing.primer_info import primers_coordinates_HXB2_outer as pcos_HXB2
from hivwholeseq.utils.mapping import get_number_reads

from hivwholeseq.sequencing.samples import load_sequencing_run, SampleSeq
from hivwholeseq.sequencing.samples import load_samples_sequenced as lss
from hivwholeseq.cluster.fork_cluster import fork_premapped_coverage as fork_self



# Functions
def plot_coverage(counts, offset_x=0, frags_pos=None, frags_pos_out=None,
                  title=None):
    '''Plot the coverage and the minor allele frequency'''
    cov = counts.sum(axis=1)
    cov_tot = cov.sum(axis=0)

    fig, ax = plt.subplots(1, 1, figsize=(11, 8))
    ax.plot(np.arange(len(cov_tot.T)) + offset_x, cov_tot.T, lw=2, c='k',
            label=read_types)
    ax.set_xlabel('Position [bases]')
    ax.set_ylabel('Coverage')
    ax.grid(True)

    # If the fragments positions are marked, plot them
    # Inner primers
    if frags_pos is not None:
        for i, frag_pos in enumerate(frags_pos.T):
            ax.plot(frag_pos + offset_x, 2 * [(0.97 - 0.03 * (i % 2)) * ax.get_ylim()[1]],
                    c=cm.jet(int(255.0 * i / len(frags_pos.T))), lw=2)

    # Outer primers
    if frags_pos_out is not None:
        for i, frag_pos in enumerate(frags_pos_out.T):
            ax.plot(frag_pos + offset_x, 2 * [(0.96 - 0.03 * (i % 2)) * ax.get_ylim()[1]],
                    c=cm.jet(int(255.0 * i / len(frags_pos_out.T))), lw=2)

    ax.set_xlim(0, 10000)

    if title is not None:
        ax.set_title(title, fontsize=18)

    plt.tight_layout(rect=(0, 0, 1, 0.95)) 


def check_premap(data_folder, adaID, fragments, seq_run, samplename,
                 qual_min=30, match_len_min=10,
                 maxreads=-1, VERBOSE=0,
                 savefig=True,
                 title=None):
    '''Check premap to reference: coverage, etc.'''
    refseq = SeqIO.read(get_reference_premap_filename(data_folder, adaID), 'fasta')

    # FIXME: do this possibly better than parsing the description!
    try:
        fields = refseq.description.split()
        refseq_start = int(fields[fields.index('(indices') - 3])
    except ValueError:
        refseq_start = 550

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
        # NOTE: In a lot of old files, it says F3o instead of F3ao
        if 'F3o' in fragtmp:
            fragtmp[fragtmp.index('F3o')] = 'F3ao'
        elif 'F3i' in fragtmp:
            fragtmp[fragtmp.index('F3i')] = 'F3ai'


        frags_pos = np.array([postmp[fragtmp.index(fr)] for fr in fragments], int).T

    else:
        frags_pos = None
    
    frags_pos_out = None

    # Open BAM and scan reads
    input_filename = get_premapped_filename(data_folder, adaID, type='bam')
    if not os.path.isfile(input_filename):
        if VERBOSE:
            print 'Premapped BAM file not found'
        return (None, None)

    # Count reads if requested
    n_reads = get_number_reads(input_filename)
    if VERBOSE:
        print 'N. of reads:', n_reads

    # Get counts
    counts, inserts = get_allele_counts_insertions_from_file_unfiltered(input_filename,
                                                             len(refseq),
                                                             qual_min=qual_min,
                                                             match_len_min=match_len_min,
                                                             maxreads=maxreads,
                                                             VERBOSE=VERBOSE)

    # Plot results
    if title is None:
        title=', '.join(['run '+seq_run+' '+adaID,
                         'sample '+samplename,
                         'reads '+str(min(maxreads, n_reads))+'/'+str(n_reads),
                        ])
    plot_coverage(counts,
                  offset_x=refseq_start,
                  frags_pos=frags_pos,
                  frags_pos_out=frags_pos_out,
                  title=title)

    if savefig:
        from hivwholeseq.sequencing.adapter_info import foldername_adapter
        plt.savefig(data_folder+foldername_adapter(adaID)+'figures/coverage_premapped_'+samplename+'.png')

    return (counts, inserts)



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Check premapped coverage',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    runs_or_samples = parser.add_mutually_exclusive_group(required=True)
    runs_or_samples.add_argument('--runs', nargs='+',
                                 help='Seq run to analyze (e.g. Tue28)')
    runs_or_samples.add_argument('--samples', nargs='+',
                                 help='Samples to analyze (e.g. 31440_PCR1')
    parser.add_argument('--adaIDs', nargs='+',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--maxreads', type=int, default=1000,
                        help='Number of reads analyzed')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--titles', nargs='*', default=None,
                        help='Give a title to the figure')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--noshow', action='store_false', dest='show',
                        help='Do not show the plot, just save it')
    parser.add_argument('--persist', action='store_true',
                        help='Plot show blocks thread')

    args = parser.parse_args()
    samplenames = args.samples
    seq_runs = args.runs
    adaIDs = args.adaIDs
    submit = args.submit
    maxreads = args.maxreads
    VERBOSE = args.verbose
    titles = args.titles
    show = args.show
    persist = args.persist

    if not show:
        plt.ioff()

    samples = lss()
    if samplenames is not None:
        samples = samples.loc[samples.index.isin(samplenames)]

    else:
        ind = np.zeros(len(samples), bool)
        for seq_run in seq_runs:
            dataset = load_sequencing_run(seq_run)
            data_folder = dataset.folder
            samples_run = dataset.samples
            # If the script is called with no adaID, iterate over all
            if adaIDs is not None:
                samples_run = samples_run.loc[samples_run.adapter.isin(adaIDs)]

            ind |= samples.index.isin(samples_run.index)

        samples = samples.loc[ind]

    for i, (samplename, sample) in enumerate(samples.iterrows()):
        sample = SampleSeq(sample)
        seq_run = sample['seq run']
        data_folder = sample.sequencing_run.folder
        adaID = sample.adapter
        fragments = sample.regions_complete

        if VERBOSE:
            print seq_run, adaID
        
        if submit:
            fork_self(samplename, maxreads=maxreads, VERBOSE=VERBOSE)
            continue

        if titles is not None:
            title = titles[i]
        else:
            title = None

        (counts, inserts) = check_premap(data_folder, adaID,
                                         fragments, seq_run, samplename,
                                         maxreads=maxreads,
                                         VERBOSE=VERBOSE,
                                         title=title)

        if show and (not submit) and (counts is not None):
            plt.ion()
            plt.show()
            if persist:
                plt.waitforbuttonpress(timeout=60)

