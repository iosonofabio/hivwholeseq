# vim: fdm=marker
'''
author:     Fabio Zanini
date:       26/10/13
content:    Study PCR-mediated recombination from the plasmid mixes.

            This script is used for:
                - PCR1 HiFi (Tue48, adaID N4-S1)
                - PCR1 Taq (Tue48, adaID N5-S1)
                - PCR2 HiFi (Tue48, adaID N6-S1)
                - PCR2 Taq (Tue48, adaID N1-S3)
'''
# Modules
import os
import sys
import argparse
import numpy as np
from collections import Counter, defaultdict
from itertools import izip
from operator import itemgetter, attrgetter
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment as MSA
import pysam
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.utils.miseq import alphal, alpha
from hivwholeseq.sequencing.samples import load_sample_sequenced
from hivwholeseq.sequencing.filenames import get_consensus_filename, get_mapped_filename, \
        get_allele_counts_filename, get_coverage_filename
from hivwholeseq.reference import load_custom_reference
from hivwholeseq.utils.mapping import align_muscle, pair_generator
from hivwholeseq.sequencing.minor_allele_frequency import filter_nus
from hivwholeseq.utils.sequence import expand_ambiguous_seq
from hivwholeseq.sequencing.filenames import get_allele_frequencies_filename
import hivwholeseq.utils.plot



# Globals
strains = ['LAI-III', '38540']



# Functions
def get_samplename(PCRtype):
    '''Get samplename from input arg'''
    samplename = 'RNA_mix_'+PCRtype[:4]+'_Taq'
    if 'Taq' not in PCRtype:
        samplename = samplename + 'HiFi'
    return samplename


def get_alignment(samples, fragment, VERBOSE=0):
    '''Load or generate alignment of mix with references'''
    import os
    from Bio import AlignIO

    names = [s.name for s in samples]
    ali_fn = samples[0].folder+'ali_'+'_'.join(names)+'_'+fragment+'.fasta'
    if os.path.isfile(ali_fn):
        ali = AlignIO.read(ali_fn, 'fasta')

    else:
        from hivwholeseq.utils.mapping import align_muscle
        sample = samples[0]
        samplerefs = samples[1:3]
        ali = align_muscle(sample.get_consensus(fragment),
                           samplerefs[0].get_consensus(fragment),
                           samplerefs[1].get_consensus(fragment),
                           sort=True)
        AlignIO.write(ali, ali_fn, 'fasta')

    return ali


def get_ends_ali(alim):
    '''Get the start/end coordinate of alignment such that all are covered'''
    has_all = (alim != '-').all(axis=0).nonzero()[0]
    return (has_all[0], has_all[-1] + 1)


def transform_coordinates_toali(coo, smat):
    '''Transform sequence into alignment coordinates'''
    coomat = (smat != '-').cumsum() - 1
    cooali = []
    for pos in coo:
        tmp = (coomat == pos).nonzero()[0]
        if len(tmp):
            cooali.append(tmp[0])
    return cooali


def transform_coordinates_fromali(coo, smat):
    '''Transform alignment into sequence coordinates'''
    coomat = (smat != '-').cumsum() - 1
    return coomat[coo]


def get_cocounts(bamfilename, ind_poly, maxreads=1000, VERBOSE=0, qual_min=30):
    '''Get the cocounts from the reads'''
    from collections import defaultdict, Counter
    from hivwholeseq.utils.mapping import pair_generator

    # The main data structure is a nested dictionary with the site pair as
    # first key, the allele pair as nested key
    counts = defaultdict(Counter)
    cocounts = defaultdict(Counter)

    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        for irp, reads in enumerate(pair_generator(bamfile)):
            if irp == maxreads:
                break

            if VERBOSE >= 2:
                if not ((irp + 1) % 1000):
                    print irp + 1

            als = {}
            for read in reads:
                pos_read = 0
                pos_ref = read.pos
                for (bt, bl) in read.cigar:
                    # Skip insertions and deletions
                    if bt == 2:
                        pos_ref += bl
                    elif bt == 1:
                        pos_read += bl

                    elif bt == 0:
                        itmp = ((ind_poly >= pos_ref) & (ind_poly < pos_ref + bl))
                        for ip_ref in ind_poly[itmp]:
                            ip_read = pos_read + (ip_ref - pos_ref)
                            if read.qual[ip_read] >= qual_min: 
                                # NOTE: we overwrite if the reads overlap. We should
                                # be safe because of high phred
                                als[ip_ref] = read.seq[ip_read]

                        pos_read += bl
                        pos_ref += bl

            # Add the alleles cocovered by this read pair to the cocounts
            pos_insert = sorted(als.keys())
            for i, pos2 in enumerate(pos_insert):
                al2 = als[pos2]
                counts[pos2][al2] += 1
                for pos1 in pos_insert[:i]:
                    al1 = als[pos1]
                    cocounts[(pos1, pos2)][(al1, al2)] += 1
    
    return counts, cocounts


def get_good_polymorphic_sites(alim, samplerefs, fragment, VERBOSE=0):
    '''Get polymorphic sites and alleles that are not gapped etc.'''
    # The mix reference might be shorter if PCR2
    (start, end) = get_ends_ali(alim)

    ind_poly = (alim != alim[0]).any(axis=0).nonzero()[0]
    ind_poly = ind_poly[(ind_poly >= start) & (ind_poly < end)]

    # Exclude gaps (we only look at SNPs, at the risk of missing some diversity)
    ind_poly = ind_poly[(alim[:, ind_poly] != '-').all(axis=0)]

    if VERBOSE >= 1:
        for i in xrange(2):
            print 'The mix consensus in fragment', fragment, 'is like', \
                  samplerefs[i].name+':', (alim[0, ind_poly] == alim[i+1, ind_poly]).sum()

    # The RNA strains are polymorphic themselves, so we need to exclude those
    # sites from the marker list
    covmin = 1000
    polymax = 1e-2
    ind_bad = set()
    for i, sampleref in enumerate(samplerefs):
        ind_bad_ref = []
        ac = sampleref.get_allele_counts(fragment, merge_read_types=True)
        for pos, ac_pos in enumerate(ac.T):
            if ac_pos.sum() < covmin:
                ind_bad_ref.append(pos)
                continue

            af_pos = 1.0 * ac_pos / ac_pos.sum()
            if np.sort(af_pos)[-2] > polymax:
                ind_bad_ref.append(pos)
                continue

        # Transform the coordinates into the alignment
        ind_bad_ali = transform_coordinates_toali(ind_bad_ref, alim[i+1])

        if VERBOSE >= 2:
            print 'Uncovered or polymorphic sites of', sampleref.name+':', ind_bad_ali, \
                  '- will exclude them'

        ind_bad |= set(ind_bad_ali)

    ind_poly = np.array(sorted(set(ind_poly) - ind_bad), int)
    
    # Get alleles and coordinates in the sequence reference
    al_poly = alim[:, ind_poly]
    ind_poly = transform_coordinates_fromali(ind_poly, alim[0])
    al_polyd = {ip: ap for (ip, ap) in izip(ind_poly, al_poly.T)}

    return al_polyd


def get_switches_conditional(cocounts, al_polyd, VERBOSE=0, freqmin=6e-3):
    '''Get the conditional switching probability'''
    switchesn = {}
    for (pos1, pos2), ald in cocounts.iteritems():
        (alp1r1, alp1r2) = al_polyd[pos1][1:]
        (alp2r1, alp2r2) = al_polyd[pos2][1:]

        ntot = (ald[(alp1r1, alp2r1)] + ald[(alp1r1, alp2r2)] +
                ald[(alp1r2, alp2r1)] + ald[(alp1r2, alp2r2)])

        if ntot < 100:
            continue

        nu11 = 1.0 * ald[(alp1r1, alp2r1)] / ntot
        nu12 = 1.0 * ald[(alp1r1, alp2r2)] / ntot
        nu21 = 1.0 * ald[(alp1r2, alp2r1)] / ntot
        nu22 = 1.0 * ald[(alp1r2, alp2r2)] / ntot

        nu1_ = nu11 + nu12
        nu_1 = nu11 + nu21
        nu2_ = nu21 + nu22
        nu_2 = nu22 + nu12

        if not (freqmin < nu1_ < 1 - freqmin):
            continue
        if not (freqmin < nu_1 < 1 - freqmin):
            continue
        
        Prandom = nu1_ * nu_2 + nu2_ * nu_1
        Pexp = nu12 + nu21
        
        switchesn[(pos1, pos2)] = Pexp / Prandom

    return switchesn



# Script
if __name__ == '__main__':

    # Parse input args: this is used to call itself recursively
    parser = argparse.ArgumentParser(description='Errors in HIV pure strains')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--maxreads', type=int, default=100,
                        help='Number of reads analyzed')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--PCRtype', default='PCR1',
                        help='Mix to study (PCR1, PCR2, PCR1Taq, PCR2Taq)')

    args = parser.parse_args()
    fragments = args.fragments
    maxreads = args.maxreads
    VERBOSE = args.verbose

    samplename = get_samplename(args.PCRtype)
    sample = load_sample_sequenced(samplename)
    adaID = sample.adapter

    samplerefs = [load_sample_sequenced(strain) for strain in strains]

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    for fragment in fragments:

        ali = get_alignment((sample, samplerefs[0], samplerefs[1]),
                            fragment,
                            VERBOSE=VERBOSE)
        alim = np.array(ali)

        al_polyd = get_good_polymorphic_sites(alim, samplerefs, fragment,
                                              VERBOSE=VERBOSE)

        if VERBOSE >= 1:
            cons = sample.get_consensus(fragment)
            print '3-way alignment in sample coordinates:'
            print '\n'.join(str(ip)+':\t'+''.join(ap)+' '+cons[ip] for ip, ap in al_polyd.iteritems())

        # Down to the reads
        # FIXME: maybe we should not use filtered ones?
        bamfilename = sample.get_mapped_filename(fragment, filtered=True)
        counts, cocounts = get_cocounts(bamfilename, np.sort(al_polyd.keys()),
                                        maxreads=maxreads,
                                        VERBOSE=VERBOSE)

        if VERBOSE >= 1:
            print 'Counts:'
            for ip, co in counts.iteritems():
                print str(ip)+':\t'+str(co)+'\t'+''.join(al_polyd[ip]),
                if al_polyd[ip][0] != co.most_common(1)[0][0]:
                    print 'Error'
                else:
                    print

        # Plot amplification bias
        if VERBOSE >= 1:
            print 'Plot amplification bias'
        fig, ax = plt.subplots()
        ax.set_xlabel('Position [bp]')
        ax.set_ylabel('Frequency of '+samplerefs[0].name)
        data = [(pos, 1.0 * co[al_polyd[pos][1]] / sum(co.itervalues()))
                for (pos, co) in counts.iteritems() if sum(co.itervalues()) > 100]
        (x, y) = np.array(data).T
        ind = x.argsort()
        x = x[ind]
        y = y[ind]
        y[y > 1-1e-4] = 1-1e-4
        y[y < 1e-4] = 1e-4
        ax.plot(x, y, lw=2, c='b')
        ax.set_title(args.PCRtype+', '+fragment)
        ax.set_ylim(0.9e-4, 1 - 0.9e-4)
        ax.set_yscale('logit')
        ax.grid(True)

        plt.tight_layout()
        plt.ion()
        plt.show()

        # Plot coverage
        if VERBOSE >= 1:
            print 'Plot coverage'
        fig, ax = plt.subplots()
        ax.set_xlabel('Distance [bp]')
        ax.set_ylabel('Coverage')
        data = [(pos, sum(co.itervalues())) for (pos, co) in counts.iteritems()]
        (x, y) = np.array(data).T
        ind = x.argsort()
        x = x[ind]
        y = y[ind]
        ax.plot(x, y, c='k', lw=2)
        ax.set_xlabel('Position [bp]')
        ax.set_ylabel('Coverage')
        ax.set_yscale('log')
        ax.set_title(args.PCRtype+', '+fragment)
        ax.grid(True)

        plt.tight_layout()
        plt.ion()
        plt.show()

        # Plot cocoverage
        if VERBOSE >= 1:
            print 'Plot cocoverage'
        fig, ax = plt.subplots()
        ax.set_xlabel('Distance [bp]')
        ax.set_ylabel('Cocoverage')
        
        keys = cocounts.keys()
        np.random.shuffle(keys)
        keys = keys[:1000]
        for (pos1, pos2) in keys:
            ald = cocounts[(pos1, pos2)]
            posm = 0.5 * (pos1 + pos2)
            color = cm.jet(1.0 * posm / alim.shape[1])
            ax.scatter(pos2 - pos1, sum(ald.itervalues()), s=30, color=color)

        ax.set_title(args.PCRtype+', '+fragment)
        ax.grid(True)

        plt.tight_layout()
        plt.ion()
        plt.show()

        # Count the switches
        if VERBOSE >= 1:
            print 'Calculate switch probability'
        switches = {}
        for (pos1, pos2), ald in cocounts.iteritems():
            alp1 = al_polyd[pos1]
            alp2 = al_polyd[pos2]
            noswitch = set([(alp1[1], alp2[1]), (alp1[2], alp2[2])])
            switch = set([(alp1[1], alp2[2]), (alp1[2], alp2[1])])
            n_noswitch = sum(ald[alpair] for alpair in noswitch)
            n_switch = sum(ald[alpair] for alpair in switch)
            n_tot = (n_switch + n_noswitch)
            if n_tot < 50:
                continue
            switches[(pos1, pos2)] = 1.0 * n_switch / n_tot

        # Plot switching probability
        if VERBOSE >= 1:
            print 'Plot switch probability'
        fig, ax = plt.subplots()
        ax.set_xlabel('Distance [bp]')
        ax.set_ylabel('Switch probability')
        
        keys = switches.keys()
        np.random.shuffle(keys)
        keys = keys[:1000]
        for (pos1, pos2) in keys:
            freq = switches[(pos1, pos2)]
            posm = 0.5 * (pos1 + pos2)
            color = cm.jet(1.0 * posm / alim.shape[1])
            ax.scatter(pos2 - pos1, freq + 1e-6, s=30, color=color)

        ax.set_title(args.PCRtype+', '+fragment)
        ax.grid(True)

        plt.tight_layout()
        plt.ion()
        plt.show()


        if VERBOSE >= 1:
            print 'Calculate switch probability, normalized'
        freqmin = 6e-3
        switchesn = get_switches_conditional(cocounts, al_polyd,
                                             VERBOSE=VERBOSE, freqmin=freqmin)


        if VERBOSE >= 1:
            print 'Plot switch probability, normalized'
        fig, ax = plt.subplots()
        ax.set_xlabel('Distance [bp]')
        ax.set_ylabel('Conditional switch probability')
        
        for (pos1, pos2), freq in switchesn.iteritems():
            posm = 0.5 * (pos1 + pos2)
            color = cm.jet(1.0 * posm / alim.shape[1])
            ax.scatter(pos2 - pos1, freq + 1e-6, s=30, color=color)

        ax.set_ylim(-0.05, 1.1)
        ax.set_title(args.PCRtype+', '+fragment)
        ax.grid(True)

        plt.tight_layout()
        plt.ion()
        plt.show()


