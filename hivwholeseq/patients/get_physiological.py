# vim: fdm=marker
'''
author:     Fabio Zanini
date:       12/06/14
content:    Plot physiological parameters, e.g. viral load, CD4+ counts.
'''
# Modules
import argparse
import numpy as np
import pandas as pd
import matplotlib.dates
import matplotlib.pyplot as plt

from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.patients.filenames import get_physiological_figure_filename



# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Check patient samples')
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--save', action='store_true',
                        help='Save the plot')

    args = parser.parse_args()
    pnames = args.patients
    VERBOSE = args.verbose
    use_save = args.save

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[patients.index.isin(pnames)]

    # Get patient and the sequenced samples
    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        if VERBOSE:
            print pname, patient.code

        samplenames = patient.samples.index
        dates = patient.dates -patient['infect date']
        times = np.array(dates).astype('float') / 24 / 3600e9
        viral_load = patient.viral_load
        counts_CD4 = patient.samples['CD4+ count']

        fig, ax1 = plt.subplots(figsize=(12, 10))
        ax2 = ax1.twinx()

        ind = np.array(viral_load > 0)
        ax1.plot(times[ind], viral_load.loc[ind], lw=2, color='steelblue')

        ind = np.array(counts_CD4 > 0)
        ax2.plot(times[ind], counts_CD4.loc[ind], lw=2, color='brown')

        ax1.set_ylim(ymin=50, ymax=viral_load.max() * 2)
        ax1.set_yscale('log')
        ax1.set_ylabel('Viral load [virions / ml plasma]')
        ax2.set_ylim(ymax=counts_CD4.max() * 1.2)
        ax2.set_ylabel('CD4+ cells [cells / ml plasma]', rotation=270)
        ax1.set_xlabel('Time [days from infection]')
        ax1.set_title(patient.code)
        ax1.grid(True)

        plt.tight_layout()

        if use_save:
            fig.savefig(get_physiological_figure_filename(pname))

    plt.ion()
    plt.show()
