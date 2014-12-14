# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/09/14
content:    Get divergence and diversity of the patient.
'''
# Modules
import os
import sys
import argparse
from operator import itemgetter
import numpy as np
from Bio import SeqIO
from matplotlib import cm
import matplotlib.pyplot as plt

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.filenames import get_allele_frequency_trajectories_filename, \
        get_allele_count_trajectories_filename
from hivwholeseq.patients.one_site_statistics import \
        plot_allele_frequency_trajectories_from_counts as plot_nus_from_act
from hivwholeseq.patients.one_site_statistics import \
        plot_allele_frequency_trajectories_from_counts_3d as plot_nus_from_act_3d
from hivwholeseq.patients.one_site_statistics import get_allele_count_trajectories
from hivwholeseq.fork_cluster import fork_get_allele_frequency_trajectory as fork_self


# Functions
def get_divergence_diversity_sliding(aft, block_length, VERBOSE=0):
    '''Get local divergence and diversity in a sliding window'''
    cons_ind = Patient.get_initial_consensus_noinsertions(aft)
    aft_nonanc = 1.0 - aft[:, cons_ind, np.arange(aft.shape[2])]
    aft_var = (aft * (1 - aft)).sum(axis=1)

    struct = np.ones(block_length)

    dg = np.ma.array(np.apply_along_axis(lambda x: np.convolve(x, struct, mode='valid'),
                                         axis=1, arr=aft_nonanc), hard_mask=True)
    ds = np.ma.array(np.apply_along_axis(lambda x: np.convolve(x, struct, mode='valid'),
                                         axis=1, arr=aft_var), hard_mask=True)

    # NOTE: normalization happens based on actual coverage
    norm = np.apply_along_axis(lambda x: np.convolve(x, struct, mode='valid'),
                               axis=1, arr=(-aft[:, 0].mask))

    dg.mask = norm < block_length
    dg /= norm

    ds.mask = norm < block_length
    ds /= norm

    x = np.arange(dg.shape[1]) + (block_length - 1) / 2.0

    return (x, dg, ds)


def get_divergence_diversity_blocks(aft, block_length, VERBOSE=0):
    '''Get local divergence and diversity in blocks'''
    cons_ind = aft[0].argmax(axis=0)
    n_blocks = aft.shape[2] // block_length
    dg = np.zeros((len(times), n_blocks))
    ds = np.zeros_like(dg)
    for n_block in xrange(n_blocks):
        for pos in xrange(block_length):
            pos += n_block * block_length
            af = aft[:, :, pos]
            af_nonanc = 1.0 - af[:, cons_ind[pos]]
            dg[:, n_block] += af_nonanc
            ds[:, n_block] += (af * (1 - af)).sum(axis=1)
    dg /= block_length
    ds /= block_length

    x = (np.arange(n_blocks) + 0.5) * block_length

    return (x, dg, ds)





# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get local divergence and diversity',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+', default=['all'],
                        help='Patients to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the result')
    parser.add_argument('--block-length', type=int, default=150,
                        help='Length of block to consider')
    parser.add_argument('--sliding', action='store_true',
                        help='Use a sliding window instead of dividing into blocks')
    parser.add_argument('--include-coverage', action='store_true', dest='include_cov',
                        help='Also calculate/plot the coverage')
    parser.add_argument('--save', action='store_true',
                        help='Save the divergence and diversity to file')

    args = parser.parse_args()
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    plot = args.plot
    block_length = args.block_length
    use_sliding = args.sliding
    use_coverage = args.include_cov
    save_to_file = args.save

    patients = load_patients()
    if pnames != ['all']:
        patients = patients.iloc[patients.index.isin(pnames)]

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 2:
        print 'fragments', fragments

    dgs = {}
    dss = {}

    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        for ifr, fragment in enumerate(fragments):
            if VERBOSE >= 1:
                print pname, fragment

            # NOTE: we should use depth instead of cov_min but I need to resync times
            aft, ind = patient.get_allele_frequency_trajectories(fragment,
                                                                 cov_min=100)

            # NOTE: Ns should be excluded from diversity and divergence
            aft = aft[:, :5, :]

            if use_coverage:
                (covt, ind2) = patient.get_coverage_trajectories(fragment)
                if set(ind).symmetric_difference(set(ind2)):
                    raise ValueError('Indices for allele freqs and coverage differ!')

            if use_sliding:
                (x, dg, ds) = get_divergence_diversity_sliding(aft, block_length,
                                                               VERBOSE=VERBOSE)
            else:
                (x, dg, ds) = get_divergence_diversity_blocks(aft, block_length,
                                                              VERBOSE=VERBOSE)

            # FIXME: avoid this var to get different conv and aft indices
            times = patient.times[ind]
            dgs[(pname, fragment)] = (patient.times[ind], dg)
            dss[(pname, fragment)] = (patient.times[ind], ds)

            if save_to_file:
                from hivwholeseq.patients.filenames import \
                        get_divergence_trajectories_local_filename, \
                        get_diversity_trajectories_local_filename

                # savez does not support masked arrays
                dg_save = dg.data.copy()
                dg_save[dg.mask] = -1
                ds_save = ds.data.copy()
                ds_save[ds.mask] = -1

                fn_out = get_divergence_trajectories_local_filename(pname, fragment)
                np.savez(fn_out, ind=ind, dg=dg_save, block_length=[block_length])
                fn_out = get_diversity_trajectories_local_filename(pname, fragment)
                np.savez(fn_out, ind=ind, ds=ds_save, block_length=[block_length])
                if VERBOSE >= 1:
                    print 'saved to file'


            if plot:
                fig, axs = plt.subplots(1, 2 + use_coverage,
                                        figsize=(15 + 4 * use_coverage, 8))
                fig.suptitle(pname+', '+fragment+\
                             ' (block length = '+str(block_length)+')', fontsize=18)
                ax1 = axs[0]
                ax1.set_xlabel('Position [bp]')
                ax1.set_ylabel('Divergence')
                ax1.set_yscale('log')
                ax1.set_ylim(0.9e-4, 5e-1)
                ax1.set_xlim(-150, aft.shape[2] + 50)
                ax1.grid(True)

                ax2 = axs[1]
                ax2.set_xlabel('Position [bp]')
                ax2.set_ylabel('Diversity')
                ax2.set_yscale('log')
                ax2.set_ylim(0.9e-4, 5e-1)
                ax2.set_xlim(-150, aft.shape[2] + 50)
                ax2.grid(True)

                if use_coverage:
                    ax3 = axs[2]
                    ax3.set_xlabel('Position [bp]')
                    ax3.set_ylabel('Coverage')
                    ax3.set_yscale('log')
                    ax3.set_ylim(1e-1, 7e5)
                    ax3.set_xlim(-150, aft.shape[2] + 50)
                    ax3.grid(True)

                for it in xrange(len(patient.times[ind])):
                    ax1.plot(x, dg[it] + 1e-4,
                             lw=2, label=str(patient.times[ind][it]),
                             color=cm.jet(1.0 * it / len(patient.times[ind])),
                             alpha=0.5)

                    ax2.plot(x, ds[it] + 1e-4,
                             lw=2,
                             color=cm.jet(1.0 * it / len(patient.times[ind])),
                             alpha=0.5)

                    if use_coverage:
                        ax3.plot(np.arange(aft.shape[2]), covt[it],
                                 lw=2, label=patient.samples[ind].index[it],
                                 color=cm.jet(1.0 * it / len(times)))

                # Plot averages
                ax1.scatter([-100] * len(times), dg.mean(axis=1), s=40,
                            c=[cm.jet(1.0 * it / len(times)) for it in xrange(len(times))],
                            edgecolor='none', alpha=0.5)

                ax2.scatter([-100] * len(times), ds.mean(axis=1), s=40,
                            c=[cm.jet(1.0 * it / len(times)) for it in xrange(len(times))],
                            edgecolor='none', alpha=0.5)

                if use_coverage:
                    ax3.scatter([-100] * len(times), covt.mean(axis=1), s=40,
                                c=[cm.jet(1.0 * it / len(times)) for it in xrange(len(times))],
                                edgecolor='none', alpha=0.5)

                ax1.legend(loc=4, fontsize=10, ncol=(1 + ((len(times) - 1) // 3)),
                           title='Time from transmission [days]')

                if use_coverage:
                    ax3.legend(loc=3, fontsize=10, ncol=(1 + ((len(times) - 1) // 3)),
                               title='Sample name:')

                if fragment == 'genomewide':
                    refann = patient.get_reference(fragment, 'gb')
                    for fea in refann.features:
                        if fea.type != 'fragment':
                            continue

                        frn = int(fea.id[1])
                        frn_even = 1 - (frn % 2)
                        x = (fea.location.nofuzzy_start, fea.location.nofuzzy_end)
                        y = [0.95**(1 + frn_even) * ax1.get_ylim()[1]] * 2
                        ax1.plot(x, y, color='k', lw=2)

                plt.tight_layout(rect=(0, 0, 1, 0.95))
                #plt.savefig('/ebio/ag-neher/home/fzanini/phd/sequencing/figures/'+\
                #            'divergence_diversity_local_'+pname+'_'+fragment+'.png')
                #plt.close()
       
    if plot:   
        plt.ion()
        plt.show()

