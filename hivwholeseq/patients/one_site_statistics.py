# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/02/14
content:    Collection of functions to do single site statistics (allele counts,
            coverage, allele frequencies) on patients.
'''
# Modules
import os
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio import SeqIO

from hivwholeseq.miseq import alpha



# Functions
def get_allele_count_trajectories(pname, samplenames, fragment, use_PCR1=1,
                                  VERBOSE=0):
    '''Get allele counts for a single patient sample'''
    if VERBOSE >= 1:
        print 'Getting allele counts:', pname, fragment

    from hivwholeseq.patients.filenames import get_initial_reference_filename, \
            get_allele_counts_filename
    from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file

    refseq = SeqIO.read(get_initial_reference_filename(pname, fragment), 'fasta')
    fns = []
    samplenames_out = []
    for samplename_pat in samplenames:

        # PCR1 filter here
        fn1 = get_allele_counts_filename(pname, samplename_pat, fragment, PCR=1)
        fn2 = get_allele_counts_filename(pname, samplename_pat, fragment, PCR=2)
        if use_PCR1 == 0:
            for PCR, fn in enumerate((fn1, fn2), 1):
                if os.path.isfile(fn):
                    fns.append(fn)
                    samplenames_out.append((samplename_pat, PCR))
                    if VERBOSE >= 3:
                        print samplename_pat, PCR

        elif use_PCR1 == 1:
            if os.path.isfile(fn1):
                fns.append(fn1)
                samplenames_out.append((samplename_pat, 1))
                if VERBOSE >= 3:
                    print samplename_pat, 1
            elif os.path.isfile(fn2):
                fns.append(fn2)
                samplenames_out.append((samplename_pat, 2))
                if VERBOSE >= 3:
                    print samplename_pat, 2
        elif use_PCR1 == 2:
            if os.path.isfile(fn1):
                fns.append(fn1)
                samplenames_out.append((samplename_pat, 1))
                if VERBOSE >= 3:
                    print samplename_pat, 1

    act = np.zeros((len(fns), len(alpha), len(refseq)), int)
    for i, fn in enumerate(fns):
        # Average directly over read types?
        act[i] = np.load(fn).sum(axis=0)

    return (samplenames_out, act)


def get_allele_frequency_trajectories(pname, samples, fragment, qual_min=30, VERBOSE=0):
    '''Scan the reads of all samples and write to a single file'''
    if VERBOSE >= 1:
        print 'Getting allele frequency trajectories:', pname, fragment

    from hivwholeseq.patients.filenames import get_initial_reference_filename, \
            get_mapped_to_initial_filename, get_allele_frequency_trajectories_filename, \
            get_allele_count_trajectories_filename
    from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file, \
            get_allele_counts_insertions_from_file_unfiltered, \
            filter_nus

    refseq = SeqIO.read(get_initial_reference_filename(pname, fragment), 'fasta')

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
                                       ntemplates=None, figax=None, offset=True):
    '''Plot the allele frequency trajectories from a patient'''
    if logit:
        import hivwholeseq.plot_utils
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import numpy as np

    if figax is not None:
        fig, ax = figax
    else:
        fig, ax = plt.subplots(1, 1)

    if offset:
        offset_step = 10.0 / nus.shape[2] # 10 days max offset
    else:
        offset_step = 0

    for i in xrange(nus.shape[2]):
        for j in xrange(nus.shape[1]):
            nu = nus[:, j, i]
            mask = nus.mask[:, j, i]
            if mask.all():
                continue
            
            times_nu = times[-mask]
            nu = nu[-mask]

            if nu[0] > 0.5:
                continue

            if (nu < threshold).all():
                continue

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

            ax.plot(times_nu + i * offset_step, nu, lw=1.5, ls=ls,
                    color=cm.jet(int(255.0 * i / nus.shape[2])))

    if ntemplates is not None:
        depthmax = 1.0 / ntemplates
        y = depthmax
        ax.plot(times, y, lw=3.5, ls='-',
                color='k', label='Max depth (# templates)')
        ax.plot(times, 1 - y, lw=3.5, ls='-',
                color='k')

    # Plot PCR error threshold
    y = np.repeat(2e-3, len(times))
    ax.plot(times, y, lw=3.5, ls='-',
            color='k', label='Max depth (PCR error rate)')
    ax.plot(times, 1 - y, lw=3.5, ls='-',
            color='k')

    ax.set_xlim(times[0] -10, times[-1] + 10)
    ax.set_xlabel('Time [days from initial sample]')
    if logit:
        ax.set_yscale('logit')
        ax.set_ylim(10**(-4.1), 1 - 10**(-4.1))
    else:
        ax.set_yscale('log')
        ax.set_ylim(9e-5, 1.5)

    ax.set_ylabel(r'$\nu$', fontsize=16)
    ax.set_title(title)
    ax.grid(True)

    return (fig, ax)


def plot_allele_frequency_trajectories_3d(times, nus, title='', VERBOSE=0,
                                          threshold=0.1, logit=False,
                                          ntemplates=None, figax=None):
    '''Plot the allele freq traj in 3D'''
    # NOTE: logscales in 3d plots is an open matplotlib issue (!!), forget logit ones
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    import matplotlib.pyplot as plt
    import numpy as np

    if figax is not None:
        fig, ax = figax
    else:
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

    if ntemplates is not None:
        # It's a poly-plane in 3d, parallel to Y axis
        X, Y = np.meshgrid([times[0] - 1] + list(times) + [times[-1] + 1],
                           [0, nus.shape[2] - 1])
        Zs = [np.tile([1.0 / ntemplates[0]] + list(1.0 / ntemplates) +
                      [1.0 / ntemplates[-1]], (2, 1))]
        Zs.append(1 - Zs[0])
        for Z in Zs:
            if logit:
                Z = np.log10((Z / (1.0 - Z)))
            else:
                Z = np.log10(Z)

            ax.plot_surface(X, Y, Z,
                            rstride=1, cstride=1, color='k',
                            linewidth=0, antialiased=False)

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

    return (fig, ax)


def plot_allele_frequency_trajectories_from_counts(times, act, title='', VERBOSE=0,
                                       threshold=0.1, options=[], logit=False,
                                       ntemplates=None, figax=None):
    '''Plot the allele frequency trajectories from a patient'''
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import numpy as np

    if figax is not None:
        fig, ax = figax
    else:
        fig, ax = plt.subplots(1, 1)

    for i in xrange(act.shape[2]):
        cov = act[:, :, i].sum(axis=1)
        # Take only time points with enough coverage
        ind = (cov > 100)
        if not ind.sum():
            continue

        for j in xrange(act.shape[1]):
            t = times[ind]
            nu = 1.0 * act[ind, j, i] / cov[ind]

            if (nu[0] < 0.5) and (nu > threshold).any():

                # Use dashed lines for synonymous if requested
                if 'syn-nonsyn' in options:
                    cod_initial = alpha[act[0, :, i - i%3: i - i%3 + 3].argmax(axis=0)]
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
                ax.plot(t, y, lw=1.5, ls=ls,
                        color=cm.jet(int(255.0 * i / act.shape[2])))

    if ntemplates is not None:
        depthmax = 1.0 / ntemplates
        if logit:
            y1 = np.log10(depthmax/(1 - depthmax))
            y2 = np.log10((1 - depthmax)/ depthmax)
            ax.plot(times, y1, lw=3.5, ls='-',
                    color='k', label='Max depth (# templates)')
            ax.plot(times, y2, lw=3.5, ls='-',
                    color='k')
        else:
            y = depthmax
            ax.plot(times, y, lw=3.5, ls='-',
                    color='k', label='Max depth (# templates)')

    ax.set_xlim(times[0] -10, times[-1] + 10)
    ax.set_xlabel('Time [days from transmission]')
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

    return (fig, ax)


def plot_allele_frequency_trajectories_from_counts_3d(times, act, title='', VERBOSE=0,
                                          threshold=0.1, options=[], logit=False,
                                                     figax=None):
    '''Plot the allele freq traj in 3D'''
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    import matplotlib.pyplot as plt
    import numpy as np

    if figax is not None:
        fig, ax = figax
    else:
        fig = plt.figure(figsize=(16, 12))
        ax = fig.gca(projection='3d')
        ax.view_init(5, 150)

    for i in xrange(act.shape[2]):
        cov = act[:, :, i].sum(axis=1)
        # Take only time points with enough coverage
        ind = (cov > 100)
        if not ind.sum():
            continue

        for j in xrange(act.shape[1]):
            t = times[ind]
            nu = 1.0 * act[ind, j, i] / cov[ind]

            if (nu[0] < 0.5) and (nu > threshold).any():

                # Use dashed lines for synonymous if requested
                if 'syn-nonsyn' in options:
                    cod_initial = alpha[act[0, :, i - i%3: i - i%3 + 3].argmax(axis=0)]
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
                    ax.plot(t, [i] * len(t), np.log10((nu + 1e-4)/(1-1e-4-nu)),
                            lw=2, ls=ls,
                            color=cm.jet(int(255.0 * i / act.shape[2])))
                else:
                    ax.plot(t, [i] * len(t), np.log10(nu + 1e-4),
                            lw=2, ls=ls,
                            color=cm.jet(int(255.0 * i / act.shape[2])))

    ax.set_xlim(times[0] -10, times[-1] + 10)
    ax.set_xlabel('Time [days from initial sample]')
    ax.set_ylabel('Position [bp]')
    ax.set_title(title)

    if logit:
        ax.set_zlim(-4.1, 4.1)
        trfun = lambda x: np.log10(x / (1 - x))
        tickloc = np.array([0.0001, 0.01, 0.5, 0.99, 0.9999])
        ax.set_zticks(trfun(tickloc))
        ax.set_zticklabels(map(str, tickloc))
        from matplotlib.ticker import FixedLocator
        ticklocminor = np.concatenate([[10**po * x for x in xrange(2 , 10)] for po in xrange(-4, -1)] + \
                                      [[0.1 * x for x in xrange(2 , 9)]] + \
                                      [[1 - 10**po * (10 - x) for x in xrange(2, 10)]
                                       for po in xrange(-2, -5, -1)])
        ax.zaxis.set_minor_locator(FixedLocator(trfun(ticklocminor)))
        ax.set_zlabel(r'$\nu$', fontsize=18)
    else:
        ax.set_zlim(-4.1, 0.1)
        ax.set_zlabel(r'$\log_{10} \nu$', fontsize=18)

    return (fig, ax)


def plot_coverage_trajectories(times, covt, title='', labels=None, legendtitle='',
                               VERBOSE=0):
    '''Plot coverage over time'''
    from matplotlib import cm
    import matplotlib.pyplot as plt
    import numpy as np

    if labels is None:
        labels = map(str, times)
        if not legendtitle:
            legendtitle = 'Time from transmission [days]:'

    fig, ax = plt.subplots(figsize=(16, 8))
    for it, t in enumerate(times):
        ax.plot(np.arange(covt.shape[1]), covt[it] + 0.1,
                lw=2,
                c=cm.jet(1.0 * it / len(times)),
                label=labels[it])

    ax.set_xlabel('Position [bp]')
    ax.set_xlim(0, covt.shape[1])
    ax.set_ylabel('Coverage', fontsize=18)
    ax.set_ylim(1, 5e5)
    ax.set_yscale('log')
    ax.legend(loc=3, title=legendtitle, fontsize=14,
              ncol=(1 + (len(times) - 1) // 4))
    ax.set_title(title)
    ax.grid(True)


def plot_coverage_trajectories_3d(times, covt, title='', VERBOSE=0):
    '''Plot coverage over time'''
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    import matplotlib.pyplot as plt
    import numpy as np

    fig = plt.figure(figsize=(16, 12))
    ax = fig.gca(projection='3d')
    ax.view_init(5, 150)
    (X, Y) = np.meshgrid(times, np.arange(covt.shape[1]))
    ax.plot_surface(X, Y, np.log10(covt.T + 0.1), lw=0, cmap=cm.jet)

    ax.set_xlim(times[0] -10, times[-1] + 10)
    ax.set_xlabel('Time [days from initial sample]')
    ax.set_ylabel('Position [bp]')
    ax.set_title(title)
    ax.set_zlim(-0.2, 5)
    ax.set_zlabel(r'$\log_{10} \nu$', fontsize=18)


def plot_coverage_trajectories_heatmap(times, covt, title='', VERBOSE=0):
    '''Plot coverage over time'''
    from matplotlib import cm
    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots(figsize=(16, 4 + 0.5 * len(times)))
    h = ax.imshow(np.log10(covt + 0.1),
                  interpolation='nearest',
                  aspect='auto',
                  vmin=0, vmax=6,
                  cmap=cm.jet)
    cbar = fig.colorbar(h, ticks=np.arange(6))
    cbar.set_label(r'$\log_{10} \nu$', fontsize=18)

    ax.set_ylabel('Time [days from initial sample]')
    ax.set_xlabel('Position [bp]')
    ax.set_ylim(-0.5, len(times) - 0.5)
    ax.set_yticks(np.arange(len(times)))
    ax.set_yticklabels(map(str, times), fontsize=16)
    ax.set_title(title)

