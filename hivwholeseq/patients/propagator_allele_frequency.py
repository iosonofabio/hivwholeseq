# vim: fdm=marker
'''
author:     Fabio Zanini
date:       19/05/14
content:    Calculate and plot the propagator of allele frequencies.
'''
# Modules
import argparse
import numpy as np
import matplotlib.pyplot as plt

from hivwholeseq.patients.patients import patients as patients_all
from hivwholeseq.patients.filenames import get_initial_consensus_filename, \
        get_allele_frequency_trajectories_filename, \
        get_allele_count_trajectories_filename



# Script
if __name__ == '__main__': 


    # Parse input args
    parser = argparse.ArgumentParser(description='Get allele frequency trajectories')
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragments to analyze (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--save', action='store_true',
                        help='Save the propagator to file')
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')
    parser.add_argument('--plot', nargs='?', default=None, const='2D',
                        help='Plot the propagator')

    args = parser.parse_args()
    pnames = args.patients
    fragments = args.fragments
    VERBOSE = args.verbose
    save_to_file = args.save
    plot = args.plot
    submit = args.submit

    if pnames is None:
        patients = [p for p in patients_all if (len(set(p.times())) >= 3)]
    else:
        patients = [p for p in patients_all if p.id in pnames]

    # Prepare output structures
    bins = np.insert(np.logspace(-2.5, 0, 10), 0, 0)
    hists = {fr: np.zeros((len(bins) - 1, len(bins) - 1), float) for fr in fragments}

    for patient in patients:
        pname = patient.id
        times = patient.times()
        samplenames = patient.samples

        # Keep PCR2 only if PCR1 is absent
        ind = np.nonzero(map(lambda x: ('PCR1' in x[1]) or ((times == times[x[0]]).sum() == 1),
                             enumerate(samplenames)))[0]
        times = times[ind]
        samplenames = [s for (i, s) in enumerate(samplenames) if i in ind]
    
        # If the script is called with no fragment, iterate over all
        if not fragments:
            fragments = ['F'+str(i) for i in xrange(1, 7)]
        if VERBOSE >= 2:
            print 'fragments', fragments
    
        # Iterate over samples and fragments
        for fragment in fragments:
    
            # Submit to the cluster self if requested (--save is assumed)
            if submit:
                fork_self(pname, fragment,
                          VERBOSE=VERBOSE)
                continue
    
            if VERBOSE >= 1:
                print pname, fragment
    
            act_filename = get_allele_count_trajectories_filename(pname, fragment)
            aft_filename = get_allele_frequency_trajectories_filename(pname, fragment)
    
            aft = np.load(aft_filename)
            act = np.load(act_filename)
    
            aft[np.isnan(aft)] = 0
            aft[(aft < 1e-5) | (aft > 1)] = 0
    
            # Collect counts
            for i in xrange(aft.shape[0] - 1):
                hist = np.histogram2d(aft[i].ravel(), aft[i + 1].ravel(), bins=[bins, bins])
                hists[fragment] += hist[0]

    # Plot
    for fragment in fragments:
        fig, ax = plt.subplots()
        z = hists[fragment][1:-1, 1:-1]
        z /= z.sum(axis=1)
        z = np.log10(z)
        im = ax.imshow(z, interpolation='nearest')
        ax.set_xlabel('Initial freq')
        ax.set_ylabel('Final freq')
        ax.set_xticklabels(map('{:1.1e}'.format, np.sqrt(bins[1:-2] * bins[2:-1])), rotation=45, fontsize=10)
        ax.set_yticklabels(map('{:1.1e}'.format, np.sqrt(bins[1:-2] * bins[2:-1])), rotation=45, fontsize=10)
        plt.colorbar(im)
        ax.set_title('Propagator for allele frequencies\nbetween consecutive time points, '+fragment+'\n[log10 P(x1 | x0)]', fontsize=16)
        plt.tight_layout()

    plt.ion()
    plt.show()
