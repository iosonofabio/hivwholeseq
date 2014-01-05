# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/10/13
content:    Description module for HIV patients.
'''
# Modules
from hivwholeseq.samples import sample_table



# Globals
patients = [\
 {'id': '20097'},
 {'id': '15363'},
 {'id': '15823'},
 {'id': '15313'},
 {'id': '15376'},
 {'id': '20529'},
 ]

# Include samples in patients
for p in patients:
    tmp = sample_table.filter_patient(p['id']).sort('date')
    p['samples'] = tuple(tmp['name'])
    p['sample_table'] = tmp.copy()



# Functions
def get_patient(pname):
    '''Get the patient from the sequences ones'''
    # This is an interface function, so we can change the actual data structure
    # (efficiency is not an issue)
    from operator import itemgetter
    names = map(itemgetter('id'), patients)
    if pname in names:
        return patients[names.index(pname)]
    else:
        raise ValueError('Patient with this ID not found')


def get_sequenced_samples(patient):
    '''Get only the sequenced samples of a patient'''
    from mapping.samples import samples as samples_seq
    from mapping.samples import date_to_integer

    # Get only sequenced
    samples = [s for s in patient['samples'] if s in samples_seq]

    # Sort samples by date (for later convenience)
    samples.sort(key=lambda s: date_to_integer(samples_seq[s]['date']))

    return samples


def get_initial_sequenced_sample(patient):
    '''Get the first sequenced sample to date (order of blood sampling)'''
    return get_sequenced_samples(patient)[0]
