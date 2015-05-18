# vim: fdm=indent
'''
author:     Fabio Zanini
date:       18/05/15
content:    Get days since infection from dates.
'''
# Globals
import pandas as pd

from hivwholeseq.patients.patients import load_patients
from hivwholeseq.patients.samples import load_samples_sequenced, itersample



# Script
if __name__ == '__main__':


    samples = load_samples_sequenced(include_empty=True)
    patients = load_patients()

    dsi = []
    for samplename, sample in itersample(samples):
        t = sample.date
        t0 = patients.loc[sample.patient]['infect date best']
        dt = (t - t0).days
        dsi.append({'sample': samplename,
                    'dsi': dt,
                   })
    dsi = pd.DataFrame(dsi).set_index('sample')

    print dsi
