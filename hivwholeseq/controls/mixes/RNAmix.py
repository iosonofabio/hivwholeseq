# vim: fdm=marker
'''
author:     Fabio Zanini
date:       04/03/15
content:    Investigate the RNA mixes after PCR1 and PCR2 for in vitro recombination.
'''
# Modules
import os
import sys
import argparse
import numpy as np
import pandas as pd
import pysam
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.reference import load_custom_reference
from hivwholeseq.sequencing.samples import load_sample_sequenced



# Globals
samplenames = {'PCR1': 'MIX1_new_PCR1',
               'PCR2': 'MIX1_new_PCR2'}
refnames = ['LAI-III', '38304']
# NOTE: F1 is not very good (some shorter PCR products)
# F2 and F3 have no indels between the references, which makes life easier
# F4 and F5 have indels
# Positions that give rise to strange outliers (e.g. polymorphisms in the RNA pops)
pos_outliers_all = {'F3': set([424, 934, 872, 425, 426, 941, 302, 271, 273, 823, 441, 847, 874, 286]),
                   }



# Functions
def get_refseqs(fragment, VERBOSE=0):
    '''Get reference sequences cut to the fragment specified'''
    if VERBOSE >= 2:
        print 'Get primers to cut'
    from hivwholeseq.sequencing.primer_info import primers_PCR
    pr_fwd, pr_rev = primers_PCR[fragment]

    if VERBOSE >= 2:
        print 'Load sequences and cut primers'
    from seqanpy import align_overlap
    refseqs = {}
    for refname in refnames:
        seq = load_custom_reference(refname)

        if VERBOSE >= 3:
            print refname+':', 'start AFTER pr_fwd'
        score, ali1, ali2 = align_overlap(seq, pr_fwd, score_gapopen=-20)
        start = len(ali2.rstrip('-'))

        if VERBOSE >= 3:
            print refname+':', 'end BEFORE pr_rev'
        score, ali1, ali2 = align_overlap(seq, pr_rev, score_gapopen=-20)
        end = len(ali2) - len(ali2.lstrip('-'))

        refseqs[refname] = seq[start: end]

    if VERBOSE >= 2:
        print 'Align them'
    from hivwholeseq.utils.sequence import align_pairwise
    ali = align_pairwise(refseqs[refnames[0]],
                         refseqs[refnames[1]],
                         score_gapopen=-20)
    refseqs['ali'] = ali    

    return refseqs


def plot_strain_differences(alim, VERBOSE=0, title=''):
    '''Plot differences between the two strains of the mix'''
    fig, ax = plt.subplots()
    ind_snp = (alim[0] != alim[1]) & (alim[0] != '-') & (alim[1] != '-')
    ind_indel = (alim[0] != alim[1]) & (-ind_snp)
    indi_snp = ind_snp.nonzero()[0]
    indi_indel = ind_indel.nonzero()[0]
    if len(indi_snp):
        ax.plot(indi_snp, 1 + np.arange(len(indi_snp)), lw=2, label='SNP')
    if len(indi_indel):
        ax.plot(indi_indel, 1 + np.arange(len(indi_indel)), lw=2, label='indel')
    if title:
        ax.set_title(title)
    ax.legend(loc='upper left')
    ax.grid(True)

    plt.ion()
    plt.show()


def examine_few_reads(sample, fragment, refseqs, VERBOSE=0, maxreads=10):
    '''Examine few reads to get an idea of what it looks like'''
    from collections import Counter
    from hivwholeseq.utils.mapping import pair_generator
    from seqanpy import align_overlap

    def make_triple_alignment(aliref, aliread, refseqs):
        '''Assume no indels between refs for now'''
        alid = {'read': aliread}
        for refname in refnames:
            seq = refseqs[refname]
            seqnew = []
            pos_ref = 0
            for nuc in aliref:
                if nuc == '-':
                    seqnew.append('-')
                else:
                    seqnew.append(seq[pos_ref])
                    pos_ref += 1
            seqnew = ''.join(seqnew)
            alid[refname] = seqnew
        return alid


    def check_read(alid, n_alleles, VERBOSE=0):
        '''Check number and position of switches and alleles'''
        readseqm = np.fromstring(alid['read'], 'S1')
        refseqm1 = np.fromstring(alid[refnames[0]], 'S1')
        refseqm2 = np.fromstring(alid[refnames[1]], 'S1')

        ind_r1 = ((readseqm == refseqm1) & (readseqm != refseqm2))
        ind_r2 = ((readseqm == refseqm2) & (readseqm != refseqm1))
        ind_snp = ind_r1 | ind_r2

        pattern = ''.join(map(str, np.array(ind_r2[ind_snp], int)))
        if VERBOSE >= 3:
            print pattern

        n_alleles[refnames[0]] += ind_r1.sum()
        n_alleles[refnames[1]] += ind_r2.sum()

        n_snp = ind_snp.sum()

        # If we have too few SNPs, skip
        if n_snp < 4:
            n_switches = None
        else:

            n_switches = 0
            indi_snp = ind_snp.nonzero()[0]
            n_snp_old = 0
            n_snp_new = 1
            snp_new = ind_r2[indi_snp[0]]
            for i in indi_snp[1:]:
                # if the new SNP is like the last one, move on
                if ind_r2[i] == snp_new:
                    n_snp_new += 1

                # else, swap pointers around
                else:
                    if (n_snp_new >= 2) and (n_snp_old >= 2):
                        n_switches += 1

                    n_snp_old = n_snp_new
                    n_snp_new = 1
                    snp_new = ind_r2[i]

            # At the end of the read, check again for a last switch
            if (n_snp_new >= 2) and (n_snp_old >= 2):
                n_switches += 1

        return {'n_snps': n_snp, 'n_switches': n_switches, 'pattern': pattern}


    bamfilename = sample.get_divided_filename(sample.convert_region(fragment))

    # FIXME: assume there are no indels (valid for F2 and F3, not F4 nor F5)
    if ('-' in refseqs['ali'][0]) or ('-' in refseqs['ali'][1]):
        raise ValueError('Indels between references!')
    ref = refseqs.values()[0]

    data = []
    with pysam.Samfile(bamfilename, 'rb') as bamfile:
        for irp, read_pair in enumerate(pair_generator(bamfile)):
            if irp == maxreads:
                break

            if VERBOSE >= 2:
                if not ((irp + 1) % 50):
                    print irp + 1

            pattern = []
            n_snps = []
            n_switches = []
            n_alleles = Counter()
            for read in read_pair:
                # The read should already be trimmed for quality
                score, aliref, aliread = align_overlap(ref, read.seq,
                                                       score_gapopen=-20)

                alid = make_triple_alignment(aliref, aliread, refseqs)

                datum = check_read(alid, n_alleles, VERBOSE=VERBOSE)

                n_snps.append(datum['n_snps'])
                n_switches.append(datum['n_switches'])
                pattern.append(datum['pattern'])

            frac_alleles = {key: 1.0 * value / sum(n_alleles.itervalues())
                            for key, value in n_alleles.iteritems()}

            datum = {'frac_alleles_'+refname: frac_alleles[refname]
                     for refname in refnames}
            datum['n_snps_r1'] = n_snps[0]
            datum['n_snps_r2'] = n_snps[1]
            datum['n_switches_r1'] = n_switches[0]
            datum['n_switches_r2'] = n_switches[1]
            datum['pattern_r1'] = pattern[0]
            datum['pattern_r2'] = pattern[1]
            datum['isize'] = np.abs(read.isize)
            data.append(datum)

    return data
                

def get_cocounts_mix(sample, fragment, alim, VERBOSE=0, maxreads=100):
    '''Get the binary cocount matrix for the RNA mix'''
    from hivwholeseq.two_site_statistics import get_coallele_counts_from_file
    from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file
    from hivwholeseq.utils.sequence import alphal

    if VERBOSE >= 2:
        print 'Getting cocounts'

    bamfilename = sample.get_divided_filename(sample.convert_region(fragment))
    L = alim.shape[1]
    cocounts = get_coallele_counts_from_file(bamfilename, L,
                                           maxreads=maxreads,
                                           VERBOSE=VERBOSE,
                                          )

    counts = get_allele_counts_insertions_from_file(bamfilename, L,
                                         maxreads=maxreads,
                                         VERBOSE=VERBOSE,
                                        )[0].sum(axis=0)

    return {'counts': counts,
            'cocounts': cocounts}


def store_cocounts(data, VERBOSE=0, overwrite=False):
    '''Store cocounts to file'''
    from hivwholeseq.filenames import root_data_folder
    fn = (root_data_folder+'specific/PCR_recombination/'+
          'RNA_mix'+data['PCR']+'_cocounts_'+fragment+'.pickle')

    present = os.path.isfile(fn)
    if present and (not overwrite):
        raise ValueError('File already present')
    else:
        if VERBOSE >= 2:
            print 'Writing to file:', fn
        import cPickle as pickle
        with open(fn, 'wb') as f:
            pickle.dump(data, f, protocol=-1)


def load_cocounts(PCR, fragment):
    '''Load cocounts from file'''
    from hivwholeseq.filenames import root_data_folder
    fn = (root_data_folder+'specific/PCR_recombination/'+
          'RNA_mix'+PCR+'_cocounts_'+fragment+'.pickle')
    
    import cPickle as pickle
    with open(fn, 'rb') as f:
        return pickle.load(f)


def get_switch_vs_distance(counts, alim, VERBOSE=0):
    '''Plot the switching probability vs distance'''
    from collections import defaultdict
    Psw = defaultdict(list)

    from hivwholeseq.utils.sequence import alphal
    inuc = np.array([map(alphal.index, alim[0]),
                     map(alphal.index, alim[1])], int)

    pos_poly = (alim[0] != alim[1]).nonzero()[0]

    # Allele frequencies at polymorphic positions
    co = counts[:, :, pos_poly, pos_poly].sum(axis=0)
    co = 1.0 * co / co.sum(axis=0)

    for ipos1, pos1 in enumerate(pos_poly):
        for ipos2, pos2 in enumerate(pos_poly[:ipos1]):
            pos = [pos1, pos2]
            d = pos1 - pos2


            P = (counts[inuc[0, pos1], inuc[1, pos2], pos1, pos2] +
                 counts[inuc[1, pos1], inuc[0, pos2], pos1, pos2])

            D = (P +
                 counts[inuc[0, pos1], inuc[0, pos2], pos1, pos2] +
                 counts[inuc[1, pos1], inuc[1, pos2], pos1, pos2])

            if D < 10:
                continue

            P = 1.0 * P / D

            # Normalize by the chances, given by the allele freqs at those sites
            N = (co[inuc[0, pos1], ipos1] * co[inuc[1, pos2], ipos2] +
                 co[inuc[1, pos1], ipos1] * co[inuc[0, pos2], ipos2])

            Psw[d].append({'prob': P,
                           'prob_normalized': P / N,
                           'pos_mean': 0.5 * (pos1 + pos2),
                           'pos1': pos1,
                           'pos2': pos2,
                          })

    Psw_mean = np.ma.masked_all(600)
    Psw_std = np.ma.masked_all(600)
    for d, val in Psw.iteritems():
        if d >= len(Psw_mean):
            continue
        Psw_mean[d] = np.mean([v['prob'] for v in val])
        Psw_std[d] = np.std([v['prob'] for v in val])

    return {'mean': Psw_mean,
            'std': Psw_std,
            'full': Psw,
           }


def filter_outliers(Psw, fragment, de_novo=False):
    '''Some positions are always crazy, find them and exclude'''
    if de_novo:
        from collections import Counter
        threshold = 20
        pos_outliers = Counter()
        for d, val in Psw['full'].iteritems():
            for v in val:
                if v['prob_normalized'] > 0.2:
                    pos_outliers[v['pos1']] += 1
                    pos_outliers[v['pos2']] += 1

        # Some are real outliers, keep the rest
        pos_outliers = set([pos for pos, c in pos_outliers.iteritems() if c > threshold])
        print pos_outliers
    else:
        pos_outliers = pos_outliers_all[fragment]

    Psw_filter = {'full': {}}
    for d, val in Psw['full'].iteritems():
        val2 = []
        for v in val:
            if (v['pos1'] in pos_outliers) or (v['pos2'] in pos_outliers):
                continue
            val2.append(v)
        if len(val2):
            Psw_filter['full'][d] = val2

    return Psw_filter


def plot_switch_vs_distance(Psws, normalized=False, title=''):
    '''Plot the probability of switching vs distance'''
    import seaborn as sns
    sns.set_style('darkgrid')
    fs = 16

    fig, axs = plt.subplots(1, 2, figsize=(9, 5), sharey=True)

    if normalized:
        key = 'prob_normalized'
    else:
        key = 'prob'

    PCRs = ['PCR1', 'PCR2']

    for iax, PCR in enumerate(PCRs):
        ax = axs[iax]
        Psw = Psws[PCR]
        for d, val in Psw['full'].iteritems():
            ax.scatter([d] * len(val),
                       [v[key] for v in val],
                       c=[cm.jet(1.0 * v['pos_mean'] / 2000) for v in val],
                       marker='s',
                       s=30,
                      )

        ax.set_xlabel('Distance [bp]', fontsize=fs)
        
        if iax == 0:
            ax.set_ylabel('Switching probability', fontsize=fs)
            ax.yaxis.set_tick_params(labelsize=fs)

        ax.set_xlim(0, 600)
        ax.set_ylim(0, 1.05)
        ax.grid(True)
        ax.xaxis.set_tick_params(labelsize=fs)

        ax.set_title(PCR, fontsize=fs)

    if title:
        fig.suptitle(title, fontsize=fs+2)

    plt.tight_layout(rect=(0, 0, 1, 0.97))

    plt.ion()
    plt.show()



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Study RNA mix and in vitro recombination')
    parser.add_argument('--fragment', choices=['F'+str(i) for i in xrange(1, 7)],
                        required=True,
                        help='Fragment to study (e.g. F1 F6)')
    parser.add_argument('--PCRs', choices=('1', '2', 'both'), default='both',
                        help='PCR to study (1, 2, or both)')
    parser.add_argument('--verbose', type=int, default=2,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=100,
                        help='Max number of read pairs to analyze')

    args = parser.parse_args()
    fragment = args.fragment
    PCRs = args.PCRs
    VERBOSE = args.verbose
    maxreads = args.maxreads

    if PCRs == 'both':
        PCRs = ['PCR1', 'PCR2']
    else:
        PCRs = ['PCR'+PCRs]

    Psws = {}
    for PCR in PCRs:
        print PCR

        sample = load_sample_sequenced(samplenames[PCR])

        refseqs = get_refseqs(sample.convert_region(fragment), VERBOSE=VERBOSE)
        # NOTE: the order is the same like refnames
        alim = np.array(refseqs['ali'])

        if VERBOSE >= 3:
            print 'Differences between the strains, '+fragment
            plot_strain_differences(alim, VERBOSE=VERBOSE,
                                    title='SNPs between HIV-1 strains, '+fragment,
                                   )

        if VERBOSE >= 1:
            print 'Get cocount matrix'
        #counts = get_cocounts_mix(sample, fragment, alim, VERBOSE=VERBOSE,
        #                          maxreads=maxreads)
        #data = {'cocounts': counts['cocounts'],
        #        'counts': counts['counts'],
        #        'alim': alim,
        #        'refnames': refnames,
        #        'PCR': PCR,
        #        'fragment': fragment,
        #        }
        #store_cocounts(data, VERBOSE=VERBOSE, overwrite=True)
        data = load_cocounts(PCR, fragment)

        if VERBOSE >= 1:
            print 'Plot switch probability vs distance'
        Psw = get_switch_vs_distance(data['cocounts'], alim, VERBOSE=VERBOSE)
        Psws[PCR] = Psw

    Psws_filter = {}
    for PCR, Psw in Psws.iteritems():
        Psw_filter = filter_outliers(Psw, fragment)
        Psws_filter[PCR] = Psw_filter

    plot_switch_vs_distance(Psws_filter, title=fragment, normalized=True)

    if False:
        

        if VERBOSE >= 1:
            print 'Examine a few reads'
        data = examine_few_reads(sample, fragment, refseqs, VERBOSE=VERBOSE,
                                 maxreads=maxreads)
        data = pd.DataFrame(data)

        if VERBOSE >= 2:
            print 'Plot fraction of '+refnames[0]+' vs '+refnames[1]

            fig, ax = plt.subplots()
            ax.hist(data.loc[:, 'frac_alleles_'+refnames[0]],
                    bins=np.linspace(0, 1, 20))
            ax.set_xlabel('Fraction of '+refnames[0])
            ax.set_title('Alleles in the reads, '+fragment+', '+PCR)
            ax.set_ylim(ymax=ax.get_ylim()[1] * 1.04)
            ax.grid(True)

            plt.ion()
            plt.show()

        #NOTE: OK, the proportions in F2/F3 for sure are 50/50 (the other ones probably
        # too, just the current algorithm does not work)


        data['n_switches_sum'] = data['n_switches_r1'] + data['n_switches_r2']
        data['mixing'] = 0.5 - (data['frac_alleles_'+refnames[0]] - 0.5).abs()

        # Get the number of recombinants over the number of chances
        # (multiple switches seem to be rare)
        n_rec = (data['n_switches_sum'] > 0).sum()
        n_tot = (data['n_switches_sum'] >= 0).sum()
        cofreq = ((data['frac_alleles_'+refnames[0]] < 0.5).mean() * 
                  (data['frac_alleles_'+refnames[0]] > 0.5).mean())
        isize = data['isize'].mean()

        # The factor 2 comes because recombination can happen 0->1 and 1->0
        rho = 1.0 * n_rec / n_tot / (cofreq * 2.0) / isize
        print 'Recombination rate per bp: {:1.1e}'.format(rho)
