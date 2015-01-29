# vim: fdm=marker
'''
author:     Fabio Zanini
date:       31/03/14
content:    Measure in vitro recombination rates using two RNA templates, one from
            subtype B and one from subtype C.
'''
# Modules
import os
import sys
import argparse
from operator import itemgetter
from collections import defaultdict
import numpy as np
import pysam
from Bio import SeqIO

from hivwholeseq.miseq import alpha
from hivwholeseq.sequencing.samples import load_sequencing_run
from hivwholeseq.sequencing.filenames import get_consensus_filename, \
        get_allele_frequencies_filename, \
        get_mapped_filename
from hivwholeseq.utils.mapping import align_muscle, pair_generator


# Globals
sd = {'mixPCR1': {'run': 'Tue48', 'adaID': 'N4-S1'},
      'mixPCR2': {'run': 'Tue48', 'adaID': 'N6-S1'},
      'mixPCR1Taq': {'run': 'Tue48', 'adaID': 'N5-S1'},
      'mixPCR2Taq': {'run': 'Tue48', 'adaID': 'N1-S3'},
      #'ref1': ,
      #'ref2': ,
     }





# Script
if __name__ == '__main__':


    # Parse input args: this is used to call itself recursively
    parser = argparse.ArgumentParser(description='Errors in HIV pure strains')
    parser.add_argument('--fragment', default='F2',
                        help='Fragment to analyze (e.g. F2)')
    parser.add_argument('--maxreads', type=int, default=100,
                        help='Number of reads analyzed')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--PCRtype', default='PCR1',
                        help='Mix to study (PCR1, PCR2, PCR1Taq, PCR2Taq)')

    args = parser.parse_args()
    fragment = args.fragment
    maxreads = args.maxreads
    VERBOSE = args.verbose
    mixname = 'mix'+args.PCRtype

    # Enrich dict
    for samplename in sd:
        dataset = load_sequencing_run(sd[samplename]['run'])
        sd[samplename]['dataset'] = dataset
        sd[samplename]['folder'] = dataset['folder']

    # Get all three consensi and make an alignment
    consB = SeqIO.read(get_consensus_filename(sd['refB']['folder'],
                                              sd['refB']['adaID'], fragment), 'fasta')
    consC = SeqIO.read(get_consensus_filename(sd['refC']['folder'],
                                              sd['refC']['adaID'], fragment), 'fasta')
    consm = SeqIO.read(get_consensus_filename(sd[mixname]['folder'],
                                              sd[mixname]['adaID'], fragment), 'fasta')
    ali_cons = align_muscle(consB, consC, consm, sort=True)

    # Find private alleles in the references
    afB = np.load(get_allele_frequencies_filename(sd['refB']['folder'],
                                                  sd['refB']['adaID'], fragment))
    afC = np.load(get_allele_frequencies_filename(sd['refC']['folder'],
                                                  sd['refC']['adaID'], fragment))

    allB = map(lambda x: alpha[np.nonzero(x)[0]].tostring(), (afB > 0.01).T)
    allC = map(lambda x: alpha[np.nonzero(x)[0]].tostring(), (afC > 0.01).T)

    # NOTE: the private allele list uses coordinates in the MIX consensus! -
    # i.e. those of the reads
    # The alignment coordinates are not used after that point!
    all_priv = []
    posB = 0
    posC = 0
    pos_mix = 0
    for pos_ali in xrange(len(ali_cons[0])):
        gapB = ali_cons[0][pos_ali] == '-'
        gapC = ali_cons[1][pos_ali] == '-'
        gapm = ali_cons[2][pos_ali] == '-'

        if not gapm:
            if gapB and (not gapC):
                all_priv.append((pos_mix, '-', allC[posC]))
                posC += 1

            elif gapC and (not gapB):
                all_priv.append((pos_mix, allB[posB], '-'))
                posB += 1

            elif (not gapB) and (not gapC):
                if not (set(allB[posB]) & set(allC[posC])):
                    all_priv.append((pos_mix, allB[posB], allC[posC]))
                posB += 1
                posC += 1

            pos_mix += 1
    
        else:
            if not gapB:
                posB += 1

            if not gapC:
                posC += 1

    all_priv = np.array(all_priv, dtype=[('pos', int), ('aB', 'S4'), ('aC', 'S4')])
    all_priv_pos_list = all_priv['pos'].tolist()

    # Down to the reads
    # N of recombination events
    n_rec = 0
    # N of events/coverage per pairv (rec, bothB, bothC, tot)
    n_rec_pairs = {}
    for i, p1 in enumerate(all_priv_pos_list):
        for p2 in all_priv_pos_list[:i]:
            n_rec_pairs[(p2, p1)] = [0, 0, 0, 0]
    # N of full C reads (divided by number of observed markers)
    n_full_C = defaultdict(int)
    n_full_B = defaultdict(int)

    # Allele counts: the alphabet is: B, C, X (neither one)
    allele_counts = np.zeros((3, len(all_priv)), int)
    bamfilename = get_mapped_filename(sd[mixname]['folder'], sd[mixname]['adaID'], fragment, filtered=True)
    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        for irp, read_pair in enumerate(pair_generator(bamfile)):
            if (VERBOSE >= 2) and not ((irp + 1) % 1000):
                print irp + 1

            if irp == maxreads:
                if VERBOSE >= 2:
                    print 'Maxreads reached:', maxreads
                break

            markers_pair = []

            for read in read_pair:
                markers = {}
                pos_ref = read.pos
                pos_read = 0
                for (bt, bl) in read.cigar:
                    if bt == 1:
                        pos_read += bl
                    elif bt == 2:
                        mark_pot = all_priv[(all_priv['pos'] >= pos_ref) & (all_priv['pos'] < pos_ref + bl)]
                        for mark in mark_pot:
                            if mark['aB'] == '-':
                                markers[mark['pos']] = 'B'
                            elif mark['aC'] == '-':
                                markers[mark['pos']] = 'C'
                            else:
                                markers[mark['pos']] ='X'

                        pos_ref += bl

                    else:
                        mark_pot = all_priv[(all_priv['pos'] >= pos_ref) & (all_priv['pos'] < pos_ref + bl)]
                        for mark in mark_pot:
                            nuc = read.seq[pos_read + mark['pos'] - pos_ref]
                            if nuc in mark['aB']:
                                markers[mark['pos']] = 'B'
                            elif nuc in mark['aC']:
                                markers[mark['pos']] = 'C'
                            else:
                                markers[mark['pos']] = 'X'

                        pos_ref += bl
                        pos_read += bl

                markers_pair.append(markers)

                # Add allele counts
                for mark in markers.iteritems():
                    if mark[1] == 'B':
                        allele_counts[0, all_priv_pos_list.index(mark[0])] += 1
                    elif mark[1] == 'C':
                        allele_counts[1, all_priv_pos_list.index(mark[0])] += 1
                    else:
                        allele_counts[2, all_priv_pos_list.index(mark[0])] += 1

            # Consensus within the read: exclude markers that are not agreed upon
            markers_pair_new = {}
            for pos in markers_pair[0]:
                if (pos not in markers_pair[1]) or (markers_pair[0][pos] == markers_pair[1][pos]):
                    markers_pair_new[pos] = markers_pair[0][pos]
            for pos in markers_pair[1]:
                if (pos not in markers_pair[0]):
                    markers_pair_new[pos] = markers_pair[1][pos]
            markers_pair_new = list(markers_pair_new.iteritems())
            markers_pair_new.sort(key=itemgetter(0))

            markers_pair = markers_pair_new
            
            markers_pairstr = ''.join(map(itemgetter(1), markers_pair))
            markers_pairstr = ''.join(map(itemgetter(1), markers_pair))
            if not (('B' in markers_pairstr) and ('C' in markers_pairstr)):

                if 'B' not in markers_pairstr:
                    n_full_C[markers_pairstr.count('C')] += 1
                    if VERBOSE >= 4:
                        print '{:>8d}'.format(irp), 'FULL C', markers_pairstr

                if 'C' not in markers_pairstr:
                    n_full_B[markers_pairstr.count('B')] += 1
                    if VERBOSE >= 5:
                        print '{:>8d}'.format(irp), 'FULL B', markers_pairstr

                continue

            markers_pair_pos = map(itemgetter(0), markers_pair)
            rec_ev = []
            sold = ''
            snew = markers_pair[0][1]
            pos_susp = 0
            for i, (pos, m) in enumerate(markers_pair[1:], 1):
                if m == 'X':
                    continue
                elif m in snew:
                    snew = snew+m
                else:
                    if (len(sold) >= 2) and (len(snew) >= 2):
                        rec_ev.append(pos_susp)
                    sold = snew
                    snew = m
                    pos_susp = 0.5 * (pos + markers_pair[i-1][0])
                
            if (len(sold) >= 2) and (len(snew) >= 2):
                rec_ev.append(pos_susp)

            if not len(rec_ev):
                continue

            n_rec += len(rec_ev)

            for i, p1 in enumerate(markers_pair_pos):
                for j, p2 in enumerate(markers_pair_pos[:i]):
                    n_rec_tmp = n_rec_pairs[(p2, p1)]
                    n_rec_tmp[-1] += 1
                    for rtmp in rec_ev:
                        if p2 < rtmp < p1:
                            n_rec_tmp[0] += 1
                    
                    if (markers_pair[i][1] == 'B') and (markers_pair[j][1] == 'B'):
                        n_rec_tmp[1] += 1
                    elif (markers_pair[i][1] == 'C') and (markers_pair[j][1] == 'C'):
                        n_rec_tmp[2] += 1

            if VERBOSE >= 3:
                print '{:>8d}'.format(irp), '{:>5d}'.format(n_rec), markers_pairstr


    if VERBOSE:
        print '# of private alleles:', len(all_priv)
        print 'allele counts'
        for row in allele_counts:
            print row[:30]

    ## Store counts
    #data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/specific/PCR_recombination/'
    #allele_counts_filename = data_folder+'RNA_'+mixname+'_allele_counts_B_C_neither_'+fragment+'.npy'
    #allele_counts.dump(allele_counts_filename)

    ## Store rough estimate for rho
    #(fractions) = 1.0 * allele_counts[:2].sum(axis=1)
    #fractions /= fractions.sum()
    #rho = 1.0 * n_rec / (irp * 500 * 2 * fractions.prod())
    #print '{:1.2e}'.format(rho), 'events per base per total PCR cycles'
    #rho_rough_filename = data_folder+'rho_rough_'+mixname+'_'+fragment+'.txt'
    #with open(rho_rough_filename, 'w') as f:
    #    f.write(str(rho)+'\n')

    ## Store accurate estimate of recombination rate
    #posm = []
    #ds = []
    #rhos = []
    #rhos_norm = []
    #for ((p1, p2), (n_ev, _, _, n_tot)) in n_rec_pairs.iteritems():
    #    if n_tot > 0:
    #        posm.append(0.5 * (p2 + p1))
    #        ds.append(p2 - p1)
    #        rhos.append(1.0 * n_ev / n_tot)

    #        ac = allele_counts[:2, (all_priv['pos'] >= p1) & (all_priv['pos'] <= p2)].sum(axis=1)
    #        ac = 1.0 * ac / ac.sum()
    #        rhos_norm.append(1.0 * n_ev / (n_tot * 2 * ac.prod()))

    #posm = np.array(posm, int)
    #ds = np.array(ds, int)
    #rhos_norm = np.array(rhos_norm)
    #i = -np.isnan(rhos_norm)
    #rho = np.dot(rhos_norm[i], ds[i]) / np.dot(ds[i], ds[i])
    #print 'rho = {:1.2e}'.format(rho), 'events per base per total PCR cycles'

    #rho_filename = data_folder+'rho_'+mixname+'_'+fragment+'.txt'
    #tmp = np.array((posm, ds, rhos, rhos_norm), dtype=float).T
    #np.savetxt(rho_filename, tmp,
    #           fmt=['%5d', '%5d', '%2.3e', '%2.3e'],
    #           header='\t'.join(['# Pos (mean) [bp]', 'distance [bp]', 'rho', 'rho norm']))

    #import matplotlib.pyplot as plt
    #from matplotlib import cm
    #plt.figure()
    #plt.scatter(ds, rhos_norm, color=cm.jet(1.0 * posm / (posm).max()), alpha=0.5)
    #plt.plot([0, 300], [0, rho * 300], c='k', lw=2, label=r'$\rho = '+'{:2.2e}'.format(rho)+'$')
    #plt.title(mixname+', '+fragment)
    #plt.xlabel('Distance [bp]')
    #plt.ylabel('rho')
    #plt.tight_layout()
    #plt.ion()
    #plt.show()
