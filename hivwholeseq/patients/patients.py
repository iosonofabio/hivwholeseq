# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/10/13
content:    Description module for HIV patients.
'''
# Modules
import numpy as np

from hivwholeseq.samples import sample_table


# Classes
class patient(object):
    '''HIV patient'''
    def __init__(self, inp):
        '''Initialize a patient from a dict, a string, or another patient'''

        # Try copy constructor, then from dict, then from string
        try:
            for attr in ['id']:
                setattr(self, attr, getattr(inp, attr))

        except AttributeError:
            try:
                for attr in ['id']:
                    setattr(self, attr, inp[attr])
    
            except TypeError:
                self.id = inp
                tmp = sample_table.filter_patient(self.id).sort(['date', 'name'])
                tmp = tmp.set_index('name', drop=False)
                del tmp['patient']
                self.samples = tuple(tmp['name'])
                self.sample_table = tmp

    
    def __repr__(self):
        return 'patient(\''+self.id+'\')'
    
    
    def __str__(self):
        return 'Patient '+self.id


    def dates(self, unique=False):
        '''Get the dates of sampling'''
        dates = self.sample_table['date'].values
        if unique:
            dates = np.unique(dates)
        return dates

    
    def times(self, unique=False, subtract_initial=True):
        '''Get the times as integers'''
        from hivwholeseq.samples import date_to_integer as dti
        times = np.array(map(dti, self.dates(unique=unique)), int)
        if subtract_initial:
            times -= times.min()
        return times        


    def print_sample_table(self):
        '''Print longitudinal sample table'''
        t = self.sample_table
        for date in self.dates(unique=True):
            rows = t.iloc[(t['date'] == date).nonzero()[0]]
            del rows['description']
            print rows


    @property
    def initial_sample(self):
        '''The initial sample used as a mapping reference'''
        return self.sample_table.iloc[0]
        




# Globals
patients = [\
 patient('20097'),
 patient('15363'),
 patient('15823'),
 patient('15313'),
 patient('15376'),
 patient('20529'),
 ]

patients_dict = {p.id: p for p in patients}



# Functions
def get_patient(pname):
    '''Get the patient from the sequences ones'''
    if pname in patients_dict:
        return patients_dict[pname]
    else:
        raise ValueError('Patient with this ID not found')


def get_initial_sequenced_sample(patient):
    '''Get the first sequenced sample to date (order of blood sampling)'''
    return get_sequenced_samples(patient)[0]
