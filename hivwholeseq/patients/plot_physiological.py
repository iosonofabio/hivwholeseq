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
    from hivwholeseq.filenames import table_filename
    samples_table = pd.read_excel(table_filename, 'Samples timeline')
    patient_table = pd.read_excel(table_filename, 'Patients', index_col=0)
    if pnames is not None:
        patients = patient_table.loc[pnames]
    else:
        patients = patient_table

    # Get patient and the sequenced samples
    for p, patient in patients.iterrows():
        if VERBOSE:
            print p

        samples_p = samples_table[samples_table['patient'] == int(p)]
        dates = samples_p['date']
        viral_load = samples_p['viral load']

        datesplot = matplotlib.dates.date2num(dates)
        fig, ax = plt.subplots(1, 1)
        ax.plot_date(datesplot, viral_load, fmt='-', lw=2)
        ax.set_yscale('log')
        ax.set_ylabel('Viral load [virions / ml plasma]')
        ax.set_xlabel('Time [yrs from 1st sample]')
        ax.set_title(p)

        plt.ion()
        plt.show()
