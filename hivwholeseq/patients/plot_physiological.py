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



# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Check patient samples')
    parser.add_argument('--patients', nargs='+',
                        help='Patients to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    pnames = args.patients
    VERBOSE = args.verbose

    # FIXME: write a function for this
    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[patients.index.isin(pnames)]

    # Get patient and the sequenced samples
    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        if VERBOSE:
            print pname

        samplenames = patient.samples.index
        dates = patient.dates
        viral_load = patient.viral_load

        datesplot = matplotlib.dates.date2num(dates)
        fig, ax = plt.subplots(1, 1)
        ax.plot_date(datesplot, viral_load, fmt='-', lw=2)
        ax.set_yscale('log')
        ax.set_ylabel('Viral load [virions / ml plasma]')
        ax.set_xlabel('Time [yrs from 1st sample]')
        ax.set_title(pname)

        plt.ion()
        plt.show()
