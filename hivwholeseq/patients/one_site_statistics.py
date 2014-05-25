# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/02/14
content:    Collection of functions to do single site statistics (allele counts,
            coverage, allele frequencies) on patients.
'''
# Modules
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio import SeqIO

from hivwholeseq.miseq import alpha



# Functions
def get_allele_frequency_trajectories(pname, samples, fragment, qual_min=30, VERBOSE=0):
    '''Scan the reads of all samples and write to a single file'''
    if VERBOSE >= 1:
        print 'Getting allele frequency trajectories:', pname

    from hivwholeseq.patients.filenames import get_initial_consensus_filename, \
            get_mapped_to_initial_filename, get_allele_frequency_trajectories_filename, \
            get_allele_count_trajectories_filename
    from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file, \
            get_allele_counts_insertions_from_file_unfiltered, \
            filter_nus

    refseq = SeqIO.read(get_initial_consensus_filename(pname, fragment), 'fasta')

    # Prepare output data structures
    cos_traj = np.zeros((len(samples), len(alpha), len(refseq)), int)
    nus_traj = np.zeros((len(samples), len(alpha), len(refseq)))
    
    for it, sample in enumerate(samples):
        if VERBOSE >= 2:
            print pname, it, sample

        input_filename = get_mapped_to_initial_filename(pname, sample, fragment, type='bam')
        (counts, inserts) = get_allele_counts_insertions_from_file_unfiltered(input_filename,
                                                                   len(refseq),
                                                                   qual_min=qual_min,
                                                                   VERBOSE=VERBOSE)
        # Take the total counts, blending in the read types
        cou = counts.sum(axis=0)
        cos_traj[it] = cou

        # Take the filtered frequencies, blending in the read types
        nu = filter_nus(counts)
        nus_traj[it] = nu

    #FIXME: test, etc.

    return (cos_traj, nus_traj)


def plot_allele_frequency_trajectories(times, nus, title='', VERBOSE=0,
                                       threshold=0.1, options=[], logit=False,
                                       ntemplates=None):
    '''Plot the allele frequency trajectories from a patient'''
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import numpy as np

    fig, ax = plt.subplots(1, 1)
    for i in xrange(nus.shape[2]):
        for j in xrange(nus.shape[1]):
            nu = nus[:, j, i]
            if (nu[0] < 0.5) and (nu > threshold).any():

                # Use dashed lines for synonymous if requested
                if 'syn-nonsyn' in options:
                    cod_initial = alpha[nus[0, :, i - i%3: i - i%3 + 3].argmax(axis=0)]
                    cod_mut = cod_initial.copy()
                    cod_mut[i%3] = alpha[j]
                    if ('-' in cod_mut) or \
                       (str(Seq(''.join(cod_initial), ambiguous_dna).translate()) != \
                        str(Seq(''.join(cod_mut), ambiguous_dna).translate())):
                        ls = '-'
                    else:
                        ls= '--'
                else:
                    ls = '-'

                if logit:
                    y = np.log10((nu + 1e-4)/(1-1e-4-nu))
                else:
                    y = nu + 1e-4
                ax.plot(times, y, lw=1.5, ls=ls,
                        color=cm.jet(int(255.0 * i / nus.shape[2])))

    if ntemplates is not None:
        depthmax = 1.0 / ntemplates
        if logit:
            y = np.log10(depthmax/(1 - depthmax))
        else:
            y = depthmax
        ax.plot(times, y, lw=3.5, ls='-',
                color='k', label='Max depth (# templates)')

    ax.set_xlim(times[0] -10, times[-1] + 10)
    ax.set_xlabel('Time [days from initial sample]')
    if logit:
        ax.set_ylim(-4.1, 4.1)
        trfun = lambda x: np.log10(x / (1 - x))
        tickloc = np.array([0.0001, 0.01, 0.5, 0.99, 0.9999])
        ax.set_yticks(trfun(tickloc))
        ax.set_yticklabels(map(str, tickloc))
        from matplotlib.ticker import FixedLocator
        ticklocminor = np.concatenate([[10**po * x for x in xrange(2 , 10)] for po in xrange(-4, -1)] + \
                                      [[0.1 * x for x in xrange(2 , 9)]] + \
                                      [[1 - 10**po * (10 - x) for x in xrange(2, 10)] for po in xrange(-2, -5, -1)])
        ax.yaxis.set_minor_locator(FixedLocator(trfun(ticklocminor)))
    else:
        ax.set_ylim(9e-5, 1.5)
        ax.set_yscale('log')

    ax.set_ylabel(r'$\nu$', fontsize=16)
    ax.set_title(title)


def plot_allele_frequency_trajectories_3d(times, nus, title='', VERBOSE=0,
                                          threshold=0.1, logit=False):
    '''Plot the allele freq traj in 3D'''
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    import matplotlib.pyplot as plt
    import numpy as np

    fig = plt.figure(figsize=(12, 12))
    ax = fig.gca(projection='3d')
    ax.view_init(5, 150)

    for i in xrange(nus.shape[2]):
        for j in xrange(nus.shape[1]):
            nu = nus[:, j, i]
            if (nu[0] < 0.5) and (nu > threshold).any():
                if logit:
                    ax.plot(times, [i] * len(times), np.log10((nu + 1e-4)/(1-1e-4-nu)),
                            lw=2,
                            color=cm.jet(int(255.0 * i / nus.shape[2])))
                else:
                    ax.plot(times, [i] * len(times), np.log10(nu + 1e-4),
                            lw=2,
                            color=cm.jet(int(255.0 * i / nus.shape[2])))

    ax.set_xlim(times[0] -10, times[-1] + 10)
    ax.set_xlabel('Time [days from initial sample]')
    ax.set_ylabel('Position [bp]')
    if logit:
        ax.set_zlim(-4.1, 4.1)
        trfun = lambda x: np.log10(x / (1 - x))
        tickloc = np.array([0.0001, 0.01, 0.5, 0.99, 0.9999])
        ax.set_zticks(trfun(tickloc))
        ax.set_zticklabels(map(str, tickloc))
        from matplotlib.ticker import FixedLocator
        ticklocminor = np.concatenate([[10**po * x for x in xrange(2 , 10)] for po in xrange(-4, -1)] + \
                                      [[0.1 * x for x in xrange(2 , 9)]] + \
                                      [[1 - 10**po * (10 - x) for x in xrange(2, 10)] for po in xrange(-2, -5, -1)])
        ax.zaxis.set_minor_locator(FixedLocator(trfun(ticklocminor)))
        ax.set_zlabel(r'$\nu$', fontsize=18)
    else:
        ax.set_zlim(-4.1, 0.1)
        ax.set_zlabel(r'$\log_{10} \nu$', fontsize=18)
    ax.set_title(title)
    ax.grid(True)
