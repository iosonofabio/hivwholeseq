# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/05/14
content:    Read Excel patient table from Johanna/Lina.
'''
# Modules
from itertools import izip
import numpy as np
import pandas as pd



# Functions
def load_spreadsheet():
    '''Load the spreadsheet from Johanna'''
    import hivwholeseq
    JOBDIR = hivwholeseq.__path__[0].rstrip('/')+'/'
    sheet_fn = JOBDIR+'patients/Kandidater.xlsx'

    sheet = pd.read_excel(sheet_fn, 'Prover')
    return sheet


def load_patient_sheet_dict():
    '''Load a dict of info for each patient from the sheet'''
    sheet = load_spreadsheet()

    # Find patient columns
    pcols = {str(n): i for (i, n) in enumerate(sheet.columns) if type(n) == int}
    pnames = pcols.keys()

    # Extract info from our seq data structures
    pinfo = {}
    for (pname, icol) in pcols.iteritems():
        if pname != '15034':
            table = sheet.iloc[1:, icol: icol + 3].dropna(how='all')
        else:
            table = sheet.iloc[1:, [icol + 1, icol + 5, icol]].dropna(how='all')

        table.rename(columns=dict(zip(table.columns, ['date', 'viral load', 'samplename'])),
                     inplace=True)

        pinfo[pname] = table

    return pinfo


def enrich_patient_sheet_dict(pname, table):
    '''Enrich patient sheet dict from sequencing info'''  
    from hivwholeseq.patients.patients import patients as patients_old
    from pandas.tslib import Timestamp
    pold = [po for po in patients_old if po.id == pname][0]
    for (samplename, date) in izip(pold.samples, pold.dates()):
        date = Timestamp(date)
        if (date == table['date']).any():
            table['samplename'].iloc[list(table['date']).index(date)] = samplename.split('_')[0]
    return table


def enrich_from_patient_sheet_dict(patient):
    '''Enrich sequencing data structures from spreadsheet'''
    from pandas.tslib import Timestamp

    pinfo = load_patient_sheet_dict()
    table = pinfo[patient.id]
    dttable = table['date']
    sntable = np.array(table['samplename'], 'S10')
    vrtable = np.array(table['viral load'], float)

    viral_load = np.ma.masked_all(len(patient.samples), float)
    for i, (date, samplename) in enumerate(izip(patient.dates(), patient.samples)):
        date = Timestamp(date)
        tmpsn = sntable == samplename.split('_')[0]
        tmpdt = np.array(dttable == date, bool)
        if tmpsn.any():
            viral_load[i] = vrtable[tmpsn][0]
        elif tmpdt.any():
            #import ipdb; ipdb.set_trace()
            viral_load[i] = vrtable[tmpdt][0]


    patient.viral_load = viral_load



# Script
if __name__ == '__main__':

    sheet = load_spreadsheet()

    from hivwholeseq.patients.patients import patients as patients_old
    for patient in patients_old:
        enrich_from_patient_sheet_dict(patient)
