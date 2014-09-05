# vim: fdm=marker
'''
author:     Fabio Zanini
date:       31/03/14
content:    Measure in vitro recombination rates with emulsion PCR from plasmids.
'''
# Modules
import os
import sys
import argparse
from operator import itemgetter
from collections import defaultdict
import numpy as np
import gzip
import pysam
from Bio import SeqIO
from Bio import AlignIO

from hivwholeseq.miseq import alpha
from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.sequencing.filenames import get_consensus_filename, \
        get_allele_frequencies_filename, \
        get_mapped_filename, \
        get_read_filenames
from hivwholeseq.mapping_utils import align_muscle, pair_generator
from hivwholeseq.sequencing.adapter_info import foldername_adapter


# Globals
sd = {'mixconv1542': {'adaID': 'TS2'},
      'mixconv1968': {'adaID': 'TS5'},
      'mixem1542': {'adaID': 'TS6'},
      'mixem1968': {'adaID': 'TS12'},
     }

ref_names = ['NL4-3', 'SF162', 'F10']




# Script
if __name__ == '__main__':


    # Parse input args: this is used to call itself recursively
    parser = argparse.ArgumentParser(description='Errors in HIV pure strains')
    parser.add_argument('--maxreads', type=int, default=1000,
                        help='Number of reads analyzed')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--PCRtype', default='1542',
                        help='Mix to study (1542, 1968)')

    args = parser.parse_args()
    maxreads = args.maxreads
    VERBOSE = args.verbose
    PCRtype = args.PCRtype
    mixconvname = 'mixconv'+PCRtype
    mixemname = 'mixem'+PCRtype
    seqrun = 'Tuen6'

    dataset = MiSeq_runs[seqrun]
    data_folder = dataset['folder']

    ## Align references to build a coordinate map
    #if VERBOSE:
    #    print 'Aligning the 3 references'
    #from hivwholeseq.reference import load_custom_reference
    #ref_recs = [load_custom_reference(refn) for refn in ref_names]
    #ali = align_muscle(*ref_recs, sort=True)

    #if VERBOSE:
    #    print 'Saving alignment'
    #out_folder = data_folder+'emPCR_results/'
    #if not os.path.isdir(out_folder):
    #    os.mkdir(out_folder)
    #AlignIO.write(ali, out_folder+'ali_references.fasta', 'fasta')

    # Get alignment of references (there are no gaps!)
    out_folder = data_folder+'emPCR_results/'
    ali_rec = AlignIO.read(out_folder+'ali_references_'+PCRtype+'.fasta', 'fasta')
    alim = np.array(ali_rec)
    ali_names = ['NL4-3', 'SF162', 'F10']

    # Get polymorphic sites
    poss_poly = (-(alim == alim[0]).all(axis=0))
    poss_polyd = (-(alim == alim[0]).all(axis=0)).nonzero()[0]
    if VERBOSE >= 3:
        print 'Polymorphic positions:', poss_polyd
        if VERBOSE >= 4:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, 1)
            ax.plot(np.arange(len(poss_poly)), poss_poly, lw=2, c='k')
            ax.set_xlabel('Pos in amplicon [bp]')
            ax.set_ylabel('Polymorphic?')

            plt.ion()
            plt.show()

    n_recs = {}
    n_purs = {}
    n_strains = {}

    # Look at the mapped reads
    for samplename in (mixconvname, mixemname):
        if VERBOSE >= 2:
            print '------------------------------------------------------------'
            if samplename == mixconvname:
                print 'CONVENTIONAL PCR'
            else:
                print 'EMULSION PCR'
            print '------------------------------------------------------------'

        n_pure = 0
        n_recombinants = 0
        n_strain = [0, 0, 0]

        adaID = sd[samplename]['adaID']
        mapped_filename = data_folder+foldername_adapter(adaID)+'mapped.sam'
        # FIXME
        mapped_filename = mapped_filename.replace('.sam', '_subset.sam')
        with pysam.Samfile(mapped_filename, 'r') as samfile:
            for irp, reads in enumerate(pair_generator(samfile)):
                if irp == maxreads:
                    break

                if not reads[0].is_proper_pair:
                    continue

                alleles_pair = []
                for read in reads:
                    pos_start = read.pos
                    pos_end = pos_start + sum(bl for (bt, bl) in read.cigar if bt in (0, 2))

                    alleles_read = []
                    pos_read = 0
                    pos_ref = read.pos
                    qual = np.fromstring(read.qual, np.int8) - 33
                    for (bt, bl) in read.cigar:
                        if bt == 1:
                            pos_read += bl
                        elif bt == 2:
                            pos_ref += bl
                            # NOTE: no gaps in ref alignment!
                        else:
                            all_block_ind = (poss_polyd >= pos_ref) & (poss_polyd < pos_ref + bl)
                            if all_block_ind.sum():
                                all_block_indi = all_block_ind.nonzero()[0]
                                for pos_i in all_block_indi:
                                    pos = poss_polyd[pos_i]
                                    posb = pos - pos_ref
                                    al = read.seq[pos_read + posb]
                                    qualp = qual[pos_read + posb]
                                    if qualp < 30:
                                        continue
                                    al_bit = np.zeros(alim.shape[0], bool)
                                    for po in xrange(alim.shape[0]):
                                        if al == alim[po, pos]:
                                            al_bit[po] = True
                                    alleles_read.append((pos, al_bit))

                            pos_ref += bl
                            pos_read += bl

                    alleles_pair.append(alleles_read)

                poss1 = map(itemgetter(0), alleles_pair[0])
                poss2 = map(itemgetter(0), alleles_pair[1])
                poss_all = sorted(set(poss1 + poss2))
                alleles_ins = []
                for pos in poss_all:
                    if (pos in poss1) and (pos not in poss2):
                        alleles_ins.append(alleles_pair[0][poss1.index(pos)])
                    elif (pos not in poss1):
                        alleles_ins.append(alleles_pair[1][poss2.index(pos)])
                    else:
                        al1 = alleles_pair[0][poss1.index(pos)][1]
                        al2 = alleles_pair[1][poss2.index(pos)][1]
                        if (al1 == al2).all():
                            alleles_ins.append((pos, al1))

                if not len(alleles_ins):
                    continue

                (poss, alls) = zip(*alleles_ins)
                alls = np.vstack(alls)
                strain = alls.sum(axis=0).argmax()
                n_strain[strain] += 1

                # Is it a recombinant?
                is_rec = False
                for i, alli in enumerate(alls.T):
                    if i == strain:
                        continue

                    # Somewhat too loose criterion ;-)
                    if (alli & (- alls.T[strain])).sum() >= 7:
                        is_rec = True
                
                if is_rec:
                    n_recombinants += 1
                # Somewhat too loose criterion ;-)
                elif len(alli) >= 15:
                    n_pure += 1

                if (VERBOSE >= 3) or ((VERBOSE >= 2) and is_rec):
                    s1 = [str(p) for (p, v) in alleles_ins]
                    s2 = ['   '.join(v.tostring().replace('\x00', ' ').replace('\x01', 'x')) for (p, v) in alleles_ins]
                    print '\n'.join(map('\t'.join, zip(s1, s2)))
                    print 'Strain:', strain + 1, '('+ali_names[strain]+')'
                    print ''
                    import time
                    time.sleep(1)

                #import ipdb; ipdb.set_trace()


        n_recs[samplename] = n_recombinants
        n_purs[samplename] = n_pure
        n_strains[samplename] = n_strain

    for key in n_recs:
        print key, 'recs:', n_recs[key], 'pure:', n_purs[key], 'ratio:', 1.0 * n_recs[key] / n_purs[key], \
                'strains:', n_strains[key]

    ## Open the reads to have a look
    ## convPCR
    #data_folder = sd[mixconvname]['folder']
    #adaID = sd[mixconvname]['adaID']
    #read_fns = get_read_filenames(data_folder, adaID, gzip=True)
    #for read_fn in read_fns:
    #    with gzip.open(read_fn, 'rb') as f:
    #        read_iter = SeqIO.parse(f, 'fastq')
    #        for i, read in enumerate(read_iter):
    #            if i < 9900:
    #                continue

    #            print read.seq

    #            if i > 10000:
    #                break

    #sys.exit()

    ## emPCR
    #data_folder = sd[mixemname]['folder']
    #adaID = sd[mixemname]['adaID']
    #read_fns = get_read_filenames(data_folder, adaID, gzip=True)
    #for read_fn in read_fns:
    #    with gzip.open(read_fn, 'rb') as f:
    #        read_iter = SeqIO.parse(f, 'fastq')
    #        for i, read in enumerate(read_iter):
    #            if i < 9000:
    #                continue

    #            print read.seq

    #            if i > 10000:
    #                break



    ### Get all three consensi and make an alignment
    ##consB = SeqIO.read(get_consensus_filename(sd['refB']['folder'], sd['refB']['adaID'], fragment), 'fasta')
    ##consC = SeqIO.read(get_consensus_filename(sd['refC']['folder'], sd['refC']['adaID'], fragment), 'fasta')
    ##consm = SeqIO.read(get_consensus_filename(sd[mixname]['folder'], sd[mixname]['adaID'], fragment), 'fasta')
    ##ali_cons = align_muscle(consB, consC, consm, sort=True)

    ### Find private alleles in the references
    ##afB = np.load(get_allele_frequencies_filename(sd['refB']['folder'], sd['refB']['adaID'], fragment))
    ##afC = np.load(get_allele_frequencies_filename(sd['refC']['folder'], sd['refC']['adaID'], fragment))

    ##allB = map(lambda x: alpha[np.nonzero(x)[0]].tostring(), (afB > 0.01).T)
    ##allC = map(lambda x: alpha[np.nonzero(x)[0]].tostring(), (afC > 0.01).T)

    ### NOTE: the private allele list uses coordinates in the MIX consensus! - i.e. those of the reads
    ### The alignment coordinates are not used after that point!
    ##all_priv = []
    ##posB = 0
    ##posC = 0
    ##pos_mix = 0
    ##for pos_ali in xrange(len(ali_cons[0])):
    ##    gapB = ali_cons[0][pos_ali] == '-'
    ##    gapC = ali_cons[1][pos_ali] == '-'
    ##    gapm = ali_cons[2][pos_ali] == '-'

    ##    if not gapm:
    ##        if gapB and (not gapC):
    ##            all_priv.append((pos_mix, '-', allC[posC]))
    ##            posC += 1

    ##        elif gapC and (not gapB):
    ##            all_priv.append((pos_mix, allB[posB], '-'))
    ##            posB += 1

    ##        elif (not gapB) and (not gapC):
    ##            if not (set(allB[posB]) & set(allC[posC])):
    ##                all_priv.append((pos_mix, allB[posB], allC[posC]))
    ##            posB += 1
    ##            posC += 1

    ##        pos_mix += 1
    ##
    ##    else:
    ##        if not gapB:
    ##            posB += 1

    ##        if not gapC:
    ##            posC += 1

    ##all_priv = np.array(all_priv, dtype=[('pos', int), ('aB', 'S4'), ('aC', 'S4')])
    ##all_priv_pos_list = all_priv['pos'].tolist()

    ### Down to the reads
    ### N of recombination events
    ##n_rec = 0
    ### N of events/coverage per pairv (rec, bothB, bothC, tot)
    ##n_rec_pairs = {}
    ##for i, p1 in enumerate(all_priv_pos_list):
    ##    for p2 in all_priv_pos_list[:i]:
    ##        n_rec_pairs[(p2, p1)] = [0, 0, 0, 0]
    ### N of full C reads (divided by number of observed markers)
    ##n_full_C = defaultdict(int)
    ##n_full_B = defaultdict(int)

    ### Allele counts: the alphabet is: B, C, X (neither one)
    ##allele_counts = np.zeros((3, len(all_priv)), int)
    ##bamfilename = get_mapped_filename(sd[mixname]['folder'], sd[mixname]['adaID'], fragment, filtered=True)
    ##with pysam.Samfile(bamfilename, 'rb') as bamfile:
    ##    for irp, read_pair in enumerate(pair_generator(bamfile)):
    ##        if (VERBOSE >= 2) and not ((irp + 1) % 1000):
    ##            print irp + 1

    ##        if irp == maxreads:
    ##            if VERBOSE >= 2:
    ##                print 'Maxreads reached:', maxreads
    ##            break

    ##        markers_pair = []

    ##        for read in read_pair:
    ##            markers = {}
    ##            pos_ref = read.pos
    ##            pos_read = 0
    ##            for (bt, bl) in read.cigar:
    ##                if bt == 1:
    ##                    pos_read += bl
    ##                elif bt == 2:
    ##                    mark_pot = all_priv[(all_priv['pos'] >= pos_ref) & (all_priv['pos'] < pos_ref + bl)]
    ##                    for mark in mark_pot:
    ##                        if mark['aB'] == '-':
    ##                            markers[mark['pos']] = 'B'
    ##                        elif mark['aC'] == '-':
    ##                            markers[mark['pos']] = 'C'
    ##                        else:
    ##                            markers[mark['pos']] ='X'

    ##                    pos_ref += bl

    ##                else:
    ##                    mark_pot = all_priv[(all_priv['pos'] >= pos_ref) & (all_priv['pos'] < pos_ref + bl)]
    ##                    for mark in mark_pot:
    ##                        nuc = read.seq[pos_read + mark['pos'] - pos_ref]
    ##                        if nuc in mark['aB']:
    ##                            markers[mark['pos']] = 'B'
    ##                        elif nuc in mark['aC']:
    ##                            markers[mark['pos']] = 'C'
    ##                        else:
    ##                            markers[mark['pos']] = 'X'

    ##                    pos_ref += bl
    ##                    pos_read += bl

    ##            markers_pair.append(markers)

    ##            # Add allele counts
    ##            for mark in markers.iteritems():
    ##                if mark[1] == 'B':
    ##                    allele_counts[0, all_priv_pos_list.index(mark[0])] += 1
    ##                elif mark[1] == 'C':
    ##                    allele_counts[1, all_priv_pos_list.index(mark[0])] += 1
    ##                else:
    ##                    allele_counts[2, all_priv_pos_list.index(mark[0])] += 1

    ##        # Consensus within the read: exclude markers that are not agreed upon
    ##        markers_pair_new = {}
    ##        for pos in markers_pair[0]:
    ##            if (pos not in markers_pair[1]) or (markers_pair[0][pos] == markers_pair[1][pos]):
    ##                markers_pair_new[pos] = markers_pair[0][pos]
    ##        for pos in markers_pair[1]:
    ##            if (pos not in markers_pair[0]):
    ##                markers_pair_new[pos] = markers_pair[1][pos]
    ##        markers_pair_new = list(markers_pair_new.iteritems())
    ##        markers_pair_new.sort(key=itemgetter(0))

    ##        markers_pair = markers_pair_new
    ##        
    ##        markers_pairstr = ''.join(map(itemgetter(1), markers_pair))
    ##        markers_pairstr = ''.join(map(itemgetter(1), markers_pair))
    ##        if not (('B' in markers_pairstr) and ('C' in markers_pairstr)):

    ##            if 'B' not in markers_pairstr:
    ##                n_full_C[markers_pairstr.count('C')] += 1
    ##                if VERBOSE >= 4:
    ##                    print '{:>8d}'.format(irp), 'FULL C', markers_pairstr

    ##            if 'C' not in markers_pairstr:
    ##                n_full_B[markers_pairstr.count('B')] += 1
    ##                if VERBOSE >= 5:
    ##                    print '{:>8d}'.format(irp), 'FULL B', markers_pairstr

    ##            continue

    ##        markers_pair_pos = map(itemgetter(0), markers_pair)
    ##        rec_ev = []
    ##        sold = ''
    ##        snew = markers_pair[0][1]
    ##        pos_susp = 0
    ##        for i, (pos, m) in enumerate(markers_pair[1:], 1):
    ##            if m == 'X':
    ##                continue
    ##            elif m in snew:
    ##                snew = snew+m
    ##            else:
    ##                if (len(sold) >= 2) and (len(snew) >= 2):
    ##                    rec_ev.append(pos_susp)
    ##                sold = snew
    ##                snew = m
    ##                pos_susp = 0.5 * (pos + markers_pair[i-1][0])
    ##            
    ##        if (len(sold) >= 2) and (len(snew) >= 2):
    ##            rec_ev.append(pos_susp)

    ##        if not len(rec_ev):
    ##            continue

    ##        n_rec += len(rec_ev)

    ##        for i, p1 in enumerate(markers_pair_pos):
    ##            for j, p2 in enumerate(markers_pair_pos[:i]):
    ##                n_rec_tmp = n_rec_pairs[(p2, p1)]
    ##                n_rec_tmp[-1] += 1
    ##                for rtmp in rec_ev:
    ##                    if p2 < rtmp < p1:
    ##                        n_rec_tmp[0] += 1
    ##                
    ##                if (markers_pair[i][1] == 'B') and (markers_pair[j][1] == 'B'):
    ##                    n_rec_tmp[1] += 1
    ##                elif (markers_pair[i][1] == 'C') and (markers_pair[j][1] == 'C'):
    ##                    n_rec_tmp[2] += 1

    ##        if VERBOSE >= 3:
    ##            print '{:>8d}'.format(irp), '{:>5d}'.format(n_rec), markers_pairstr


    ##if VERBOSE:
    ##    print '# of private alleles:', len(all_priv)
    ##    print 'allele counts'
    ##    for row in allele_counts:
    ##        print row[:30]

    #### Store counts
    ###data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/specific/PCR_recombination/'
    ###allele_counts_filename = data_folder+'RNA_'+mixname+'_allele_counts_B_C_neither_'+fragment+'.npy'
    ###allele_counts.dump(allele_counts_filename)

    #### Store rough estimate for rho
    ###(fractions) = 1.0 * allele_counts[:2].sum(axis=1)
    ###fractions /= fractions.sum()
    ###rho = 1.0 * n_rec / (irp * 500 * 2 * fractions.prod())
    ###print '{:1.2e}'.format(rho), 'events per base per total PCR cycles'
    ###rho_rough_filename = data_folder+'rho_rough_'+mixname+'_'+fragment+'.txt'
    ###with open(rho_rough_filename, 'w') as f:
    ###    f.write(str(rho)+'\n')

    #### Store accurate estimate of recombination rate
    ###posm = []
    ###ds = []
    ###rhos = []
    ###rhos_norm = []
    ###for ((p1, p2), (n_ev, _, _, n_tot)) in n_rec_pairs.iteritems():
    ###    if n_tot > 0:
    ###        posm.append(0.5 * (p2 + p1))
    ###        ds.append(p2 - p1)
    ###        rhos.append(1.0 * n_ev / n_tot)

    ###        ac = allele_counts[:2, (all_priv['pos'] >= p1) & (all_priv['pos'] <= p2)].sum(axis=1)
    ###        ac = 1.0 * ac / ac.sum()
    ###        rhos_norm.append(1.0 * n_ev / (n_tot * 2 * ac.prod()))

    ###posm = np.array(posm, int)
    ###ds = np.array(ds, int)
    ###rhos_norm = np.array(rhos_norm)
    ###i = -np.isnan(rhos_norm)
    ###rho = np.dot(rhos_norm[i], ds[i]) / np.dot(ds[i], ds[i])
    ###print 'rho = {:1.2e}'.format(rho), 'events per base per total PCR cycles'

    ###rho_filename = data_folder+'rho_'+mixname+'_'+fragment+'.txt'
    ###tmp = np.array((posm, ds, rhos, rhos_norm), dtype=float).T
    ###np.savetxt(rho_filename, tmp,
    ###           fmt=['%5d', '%5d', '%2.3e', '%2.3e'],
    ###           header='\t'.join(['# Pos (mean) [bp]', 'distance [bp]', 'rho', 'rho norm']))

    ###import matplotlib.pyplot as plt
    ###from matplotlib import cm
    ###plt.figure()
    ###plt.scatter(ds, rhos_norm, color=cm.jet(1.0 * posm / (posm).max()), alpha=0.5)
    ###plt.plot([0, 300], [0, rho * 300], c='k', lw=2, label=r'$\rho = '+'{:2.2e}'.format(rho)+'$')
    ###plt.title(mixname+', '+fragment)
    ###plt.xlabel('Distance [bp]')
    ###plt.ylabel('rho')
    ###plt.tight_layout()
    ###plt.ion()
    ###plt.show()
