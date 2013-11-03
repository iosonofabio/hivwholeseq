# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/10/13
content:    Description module for HIV patients.
'''
# Globals
pat_20097 = {'id': '20097',
             'samples': ('VL96-15555',
                         'VL98-1253',
                         'VK99-2133',
                         'VK99-4204',
                         'VK00-0119',
                         'VK00-1524',
                         'VK01-2965',
                         'VK02-4452',
                         'VK03-3214',
                         'VK03-4298',
                         'VK04-3106',
                         'VK04-4187'),
             'comments': '',
            }

pat_15363 = {'id': '15363',
             'samples': ('VK07-0259',
                         'VK08-1001'),
             'comments': '',
            }

pat_15823 = {'id': '15823',
             'samples': ('04HR-1501',
                         '05HR-0269',
                         '06HR-0145',
                         '07HR-0248',
                         '08HR-0235'),
             'comments': '',
            }

pat_15313 = {'id': '15313',
             'samples': ('VK07-4778',
                         'VK08-2987',
                         'VK09-1685'),
             'comments': '',
            }

pat_15376 = {'id': '15376',
             'samples': ('VK06-6001',
                         'VK07-8262',
                         'VK08-8014'),
             'comments': '',
            }

pat_20529 = {'id': '20529',
             'samples': ('VK02-4864',
                         'VK03-0342',
                         'VK05-2685',
                         'VK06-1885',
                         'VK07-4218',
                         'VK08-6634',
                         'VK09-7738'),
             'comments': '',
            }

# List of patients
patients = (pat_20097,
            pat_15363,
            pat_15823,
            pat_15313,
            pat_15376,
            pat_20529)



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


# NOTE: some functions look a lot like instance methods ;-)
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
