# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/01/15
content:    Some samples have suspicious minor allele counts in the first 100 bases
            of F5, try to find out whether it's mismapping or what.
'''
# Modules
import os, sys
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.samples import load_sample_sequenced as lssp
from hivwholeseq.sequencing.samples import load_sample_sequenced as lss



# Functions
def scan_reads(bamfilename, ind, icons, consm):
    '''Scan the reads and find out the source of the minor counts'''
    from hivwholeseq.miseq import alphal
    import pysam

    ind = set(ind)

    ind_reads = []
    read_names = []
    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        for read in bamfile:
            if read.pos > 100:
                continue

            pos_ref = read.pos
            pos_read = 0
            for (bt, bl) in read.cigar:
                if bt == 1:
                    pos_read += bl
                elif bt ==2:
                    pos_ref += bl
                else:
                    seqb = np.fromstring(read.seq[pos_read: pos_read + bl], 'S1')
                    qualb = np.fromstring(read.qual[pos_read: pos_read + bl], np.int8) - 33
                    iab = map(alphal.index, seqb)
                    posb = np.arange(pos_ref, pos_ref + bl)
                    posrb = np.arange(pos_read, pos_read + bl)

                    indb = set(zip(*(posb, iab)))

                    indov = indb & ind

                    if indov:
                        print indov
                        ind_reads.append(indov)

                        if len(indov) > 5:
                            if read.is_reverse:
                                dlabel = 'rev'
                            else:
                                dlabel = 'fwd'
                            read_label = (read.qname, dlabel, read.pos, read.cigar)
                            read_names.append(read_label)

                    pos_read += bl
                    pos_ref += bl

    return ind_reads, read_names


def find_reads(bamfilename, read_names):
    '''Find a bunch of reads in a file'''
    import pysam

    qnames, directions, positions, cigars = zip(*read_names)
    dird = {True: 'rev', False: 'fwd'}

    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        for read in bamfile:
            if read.qname in qnames:
                itmp = qnames.index(read.qname)
                dtmp = directions[itmp]
                ptmp = positions[itmp]
                ctmp = cigars[itmp]
                if dird[read.is_reverse] == directions[itmp]:
                    print 'Before:', read.pos, read.cigar
                    print 'After:', ptmp, ctmp




# Script
if __name__ == '__main__':

    fragment = 'F5'

    samplename = 'VK03-3214'

    sample = lssp(samplename)
    cou = sample.get_allele_counts(fragment)

    x = np.tile(np.arange(cou.shape[1]), (cou.shape[0], 1))
    color = np.tile(np.arange(cou.shape[0]), (cou.shape[1], 1)).T

    fig, ax = plt.subplots(figsize=(12, 6))
    
    ax.scatter(x, cou + 0.1, lw=2, c=color)
    ax.set_xlabel('Position [bp]')
    ax.set_ylabel('Coverage')
    ax.set_xlim(-1, cou.shape[-1])
    ax.set_ylim(ymin=0.09)
    ax.set_yscale('log')
    ax.grid(True)
    ax.set_title(samplename)

    ax.set_xlim(-1, 130)

    cous = cou[:, :120]
    ind = zip(*(((cous > 4) & (cous < 100)).T).nonzero())

    icons = cou.argmax(axis=0)
    consm = alpha[icons]

    # Get the names of the suspicious reads
    bamfilename = sample.get_mapped_filtered_filename(fragment, decontaminated=True)
    ind_reads, read_names = scan_reads(bamfilename, ind, icons, consm)

    # Go back the pipeline to find the bug: undecontaminated
    bamfilename = sample.get_mapped_filtered_filename(fragment, decontaminated=False)
    find_reads(bamfilename, read_names)

    # Go back to the reads mapped to patient reference before filtering
    samplenameseq = samplename+'_PCR1'
    sampleseq = lss(samplenameseq)
    bamfilename =  sampleseq.get_mapped_to_initial_filename(fragment)
    find_reads(bamfilename, read_names)

    plt.ion()
    plt.show()
