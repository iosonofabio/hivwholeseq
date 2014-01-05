# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/10/13
content:    Description module for HIV samples from patients.

            Note: it is possible to read the Excel XLSX table with xlrd. It also
            works ok, but it's painful and error-prone code, and it assumes that
            table is up-to-date. Better write down the machine-relevant info
            manually down here (until otherwise specified).
'''
# Modules
import pandas as pd



# Classes
class SampleTable(pd.DataFrame):
    def filter_seq_run(self, run):
        '''Get only the samples from one sequencing run'''
        return self[self.run == run]


    def filter_patient(self, patient_id):
        '''Get only the samples from a specific patient'''
        return self[self.patient == patient_id]



# Globals
# Date is YYYY-M-D
sample_list = [\
 # Test dataset for the algorithms
 {'name': 'VK04-3106-test', 'run': 'test_tiny', 'adaID': 'TS2',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i')},
 {'name': '08HR-0235-test', 'run': 'test_tiny', 'adaID': 'TS4',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i')},
 {'name': 'VK07-4778-test', 'run': 'test_tiny', 'adaID': 'TS5',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i')},
 {'name': 'VK03-4298-test', 'run': 'test_tiny', 'adaID': 'TS6',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i')},
 {'name': 'VK09-7738-test', 'run': 'test_tiny', 'adaID': 'TS7',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i')},
 {'name': 'VK08-8014-test', 'run': 'test_tiny', 'adaID': 'TS12',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i')},
 {'name': 'VK08-6634-test', 'run': 'test_tiny', 'adaID': 'TS13',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i')},
 {'name': 'VK07-8262-test', 'run': 'test_tiny', 'adaID': 'TS14',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i')},
 {'name': 'VK08-2987-test', 'run': 'test_tiny', 'adaID': 'TS15',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i')},

 # Test MiSeq run
 {'name': 'NL4-3', 'run': 'Tue28', 'adaID': 'TS2',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i')},
 {'name': 'SF162', 'run': 'Tue28', 'adaID': 'TS4',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i')},
 {'name': 'F10', 'run': 'Tue28', 'adaID': 'TS7',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i')},
 {'name': '37024', 'run': 'Tue28', 'adaID': 'TS16',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i'),
  'description': 'patient 37024'},
 {'name': 'MIX1', 'run': 'Tue28', 'adaID': 'TS18',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i'),
  'description': 'MIX1 SF162 50% + NL4-3 50%'},
 {'name': 'MIX2', 'run': 'Tue28', 'adaID': 'TS19',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i'),
  'description': 'MIX2 SF162 95% + NL4-3 4.5% + F10 0.5%'},

 # Nextera test in Sweden
 # WE DO NOT KNOW THE ADAPTER IDS OF THESE TWO
 {'name': 'Nextera_HIV-8262-1', 'run': 'testnextera_Lina', 'adaID': '01',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
  'description': 'HIV-8262-1: 0.2 ng/ul'},
 {'name': 'Nextera_HIV-8262-2', 'run': 'testnextera_Lina', 'adaID': '02',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
  'description': 'HIV-8262-2: 0.033 ng/ul'},

 # First patient samples (nested PCR, with in vitro recombination)
 {'name': 'VK04-3106', 'run': 'Tue37', 'adaID': 'TS2', 'date': '2004-07-09',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i'),
  'patient': '20097'},
 {'name': '08HR-0235', 'run': 'Tue37', 'adaID': 'TS4', 'date': '2008-02-07',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i'),
  'patient': '15823'},
 {'name': 'VK07-4778', 'run': 'Tue37', 'adaID': 'TS5', 'date': '2007-07-12',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i'),
  'patient': '15313'},
 {'name': 'VK03-4298', 'run': 'Tue37', 'adaID': 'TS6', 'date': '2003-09-30',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i'),
  'patient': '20097'},
 {'name': 'VK09-7738', 'run': 'Tue37', 'adaID': 'TS7', 'date': '2009-09-24',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
  'patient': '20529'},
 {'name': 'VK08-8014', 'run': 'Tue37', 'adaID': 'TS12', 'date': '2008-10-28',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
  'patient': '15376'},
 {'name': 'VK08-6634', 'run': 'Tue37', 'adaID': 'TS13', 'date': '2008-09-11',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
  'patient': '20529'},
 {'name': 'VK07-8262', 'run': 'Tue37', 'adaID': 'TS14', 'date': '2007-11-27',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
  'patient': '15376'},
 {'name': 'VK08-2987', 'run': 'Tue37', 'adaID': 'TS15', 'date': '2008-04-21',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
  'patient': '15313'},

 ]

# Sample table using Pandas
sample_table = SampleTable(sample_list)

# Make a dictionary (deprecated)
samples = {}
for line in sample_list:
    dic = dict(line)
    name = dic.pop('name')
    samples[name] = dic



# Functions
def date_to_integer(date):
    '''Convert a date in the format YYYY-M-D into an integer'''
    import datetime as dt
    return dt.date.toordinal(dt.datetime(*map(int, date.split('-'))))
