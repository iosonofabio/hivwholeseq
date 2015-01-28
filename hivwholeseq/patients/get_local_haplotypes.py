# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/11/14
content:    Get local haplotypes from single read pairs, including insertions
            and deletions. This includes aggressive clustering to keep the
            multiple sequence alignments efficient.
'''
# Modules
import os
import argparse
from operator import itemgetter
from itertools import izip
import numpy as np
from Bio import SeqIO, AlignIO

from hivwholeseq.patients.patients import load_patient
from hivwholeseq.utils.argparse import RoiAction
from hivwholeseq.utils.sequence import build_msa_haplotypes as build_msa



# Functions
def store_alignments(alis, filename):
    '''Save alignments to file'''
    file_formats = {'stk': 'stockholm',
                    'fasta': 'fasta',
                    'phy': 'phylip-relaxed'}

    foldername = os.path.dirname(filename)
    if not os.path.isdir(foldername):
        raise IOError('Destination folder for file save not found')

    if os.path.isfile(filename):
        raise IOError('Destination file already exists on file system')

    file_format = filename.split('.')[-1]
    if file_format in file_formats:
        file_format = file_formats[file_format]
        AlignIO.write(alis, filename, file_format)

    elif file_format == 'zip':
        import StringIO
        import zipfile, zlib
        with zipfile.ZipFile(filename, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
            for i, ali in enumerate(alis):
                f = StringIO.StringIO()
                AlignIO.write(ali, f, 'fasta')
                zf.writestr(str(i+1)+'.fasta', f.getvalue())

    else:
        raise ValueError('File format not recognized')


def cluster_haplotypes(haplo, VERBOSE=0, min_abundance=1):
    '''Cluster haplotypes (trivial for now)'''
    from collections import Counter
    haploc = Counter()
    for (seq, count) in haplo.iteritems():
        if count >= min_abundance:
            haploc[seq] = count

    return haploc


def merge_read_pair(seq1, seq2):
    '''Merge two reads of a pair, assuming the second starts later'''
    from seqanpy import align_ladder
    (score, ali1, ali2) = align_ladder(seq1, seq2, score_gapopen=-20)
    end1 = len(ali1.rstrip('-'))
    start2 = len(ali2) - len(ali2.lstrip('-'))
    overlap_ali = np.vstack([np.fromstring(a[start2: end1], 'S1')
                             for a in (ali1, ali2)])

    overlap = overlap_ali[0]
    overlap[overlap_ali[0] != overlap_ali[1]] = 'N'
    overlap = overlap.tostring()

    seq = ali1[:start2] + overlap + ali2[end1:]
    return seq


def trim_read_roi(read, start, end):
    '''Trim a single read to a region of interest'''
    seq = []
    pos_ref = read.pos
    pos_read = 0
    for (bt, bl) in read.cigar:
        # Insertions cannot end a block, and are always copied verbatim
        # We implicitely accept leading-edge insertions
        if bt == 1:
            if pos_ref >= start:
                seq.append(read.seq[pos_read: pos_read + bl])
            pos_read += bl

        # Deletions can both start and end a block, but there's nothing to copy
        elif bt == 2:
            if pos_ref + bl >= end:
                break
            pos_ref += bl

        # Matches can both start and end blocks, and need copying
        elif bt == 0:
            if pos_ref + bl > start:
                start_inblock = max(0, start - pos_ref)
                end_inblock = min(bl, end - pos_ref)
                seq.append(read.seq[pos_read + start_inblock:
                                    pos_read + end_inblock])

                # We ignore trailing-edge insertions
                if pos_ref + bl >= end:
                    break

            pos_ref += bl
            pos_read += bl

    seq = ''.join(seq)
    return seq


def get_local_haplotypes(bamfilename, start, end, VERBOSE=0, maxreads=-1):
    '''Extract reads fully covering the region, discarding insertions'''
    import sys
    import pysam
    from hivwholeseq.utils.mapping import pair_generator
    from hivwholeseq.utils.mapping import extract_mapped_reads_subsample_open

    from collections import Counter
    haplotypes = Counter()

    with pysam.Samfile(bamfilename, 'rb') as bamfile:

        if maxreads == -1:
            reads_iter = pair_generator(bamfile)
        else:
            reads_iter =  extract_mapped_reads_subsample_open(bamfile, maxreads,
                                                              VERBOSE=VERBOSE,
                                                              pairs=True)

        for irp, reads in enumerate(reads_iter):
            if VERBOSE >= 2:
                if not ((irp + 1) % 10000):
                    if irp + 1 != 10000:
                        sys.stdout.write("\x1b[1A\n")
                    sys.stdout.write(str(irp + 1))
                    sys.stdout.flush()

            # Sort fwd read first: this is important because with our insert
            # size we know the fwd read starts <= the rev read
            is_fwd = reads[0].is_reverse
            reads = [reads[is_fwd], reads[not is_fwd]]

            # Check for coverage of the region
            start_fwd = reads[0].pos
            end_fwd = start_fwd + sum(bl for (bt, bl) in reads[0].cigar if bt in (0, 2))
            start_rev = reads[1].pos
            end_rev = start_rev + sum(bl for (bt, bl) in reads[1].cigar if bt in (0, 2))
            overlap_len = max(0, end_fwd - start_rev)

            # Various scenarios possible
            if start_fwd > start:
                continue

            if end_rev < end:
                continue

            # No single read covers the whole region AND (the insert has a whole
            # OR a very short overlap)
            if (end_fwd < end) and (start_rev > start) and (overlap_len < 20):
                continue

            # Now the good cases
            if (start_fwd <= start) and (end_fwd >= end):
                seq = trim_read_roi(reads[0], start, end)

            elif (start_rev <= start) and (end_rev >= end):
                seq = trim_read_roi(reads[1], start, end)

            else:
                seqs = [trim_read_roi(read, start, end) for read in reads]
                seq = merge_read_pair(*seqs)

            haplotypes[seq] += 1
            if VERBOSE >= 4:
                import ipdb; ipdb.set_trace()

    return haplotypes


def plot_haplotype_frequencies(times, hft, figax=None, title='',
                               picker=None):
    '''Plot haplotype frequencies'''
    import hivwholeseq.utils.plot
    from matplotlib import cm
    import matplotlib.pyplot as plt

    if figax is None:
        fig, ax = plt.subplots(figsize=(12, 7))
    else:
        fig, ax = figax

    # TODO: The hard part is finding an ordering
    hft_cum = hft.cumsum(axis=1)

    # Randomize colors to make them more visible
    colors = cm.jet(1.0 * np.arange(hft.shape[1]) / hft.shape[1])
    np.random.shuffle(colors)

    # Use fake zero/one for logit plots
    freqmin = 1e-6

    # Plot first line
    ax.fill_between(times, hft_cum[:, 0], freqmin + np.zeros(hft.shape[0]), color=colors[0],
                    label=str(0),
                    picker=picker)
    for i in xrange(1, hft.shape[1]):
        ax.fill_between(times, hft_cum[:, i],
                        np.minimum(1-freqmin, hft_cum[:, i - 1]), color=colors[i],
                        label=str(i),
                        picker=picker)

    ax.set_xlabel('Time from infection [days]')
    ax.set_ylabel('Haplotype frequency')
    ax.set_ylim(1e-4, 1 - 1e-4)
    ax.set_xlim(times[0], times[-1])

    if title:
        ax.set_title(title)

    return (fig, ax)
    


# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get local haplotypes',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--roi', required=True, action=RoiAction,
                        help='Region of interest (e.g. F1 300 350)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of reads analyzed per sample')
    parser.add_argument('--plot', action='store_true',
                        help='Plot local haplotype trajectories')
    parser.add_argument('--save', default=None,
                        help='Save to this filename')
    parser.add_argument('--countmin', type=int, default=1,
                        help='Minimal number of observations to keep the haplotype')

    args = parser.parse_args()
    pname = args.patient
    roi = args.roi
    VERBOSE = args.verbose
    maxreads = args.maxreads
    use_plot = args.plot
    save_path = args.save
    countmin = args.countmin

    patient = load_patient(pname)
    patient.discard_nonsequenced_samples()

    if VERBOSE >= 1:
        print patient.name, roi

    ind = []
    haplocs = []
    alis = [] 
    for it, (t, sample) in enumerate(izip(patient.times, patient.itersamples())):
        if VERBOSE >= 1:
            print t, sample.name

        if VERBOSE >= 2:
            print 'Get local haplotypes'
        try:
            haplo = sample.get_local_haplotypes(roi,
                                                VERBOSE=VERBOSE,
                                                maxreads=maxreads)
        except IOError:
            continue

        if VERBOSE >= 2:
            print 'Cluster haplotypes'
        haploc = cluster_haplotypes(haplo, VERBOSE=VERBOSE, min_abundance=countmin)

        if VERBOSE >= 2:
            print 'Build MSA'
        msa = build_msa(haploc, VERBOSE=VERBOSE,
                        label=patient.code+'_'+str(t)+'_days_')

        ind.append(it)
        haplocs.append(haploc)
        alis.append(msa)

    if save_path is not None:
        if VERBOSE >= 1:
            print 'Save to file:', save_path
        store_alignments(alis, save_path)

    if use_plot:
        import matplotlib.pyplot as plt

        if VERBOSE >= 1:
            print 'Plot'
        seqs_set = set()
        for haploc in haplocs:
            seqs_set |= set(haploc.keys())
        seqs_set = list(seqs_set)

        hct = np.zeros((len(seqs_set), len(haplocs)), int)
        for i, haploc in enumerate(haplocs):
            for seq, count in haploc.iteritems():
                hct[seqs_set.index(seq), i] = count
        hct = hct.T

        hft = (1.0 * hct.T / hct.sum(axis=1)).T

        plot_haplotype_frequencies(times, hft,
                                   title=patient.name+', '+' '.join(map(str, roi)))

        plt.ion()
        plt.show()
        
