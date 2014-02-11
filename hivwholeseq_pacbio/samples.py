# vim: fdm=marker
'''
author:     Fabio Zanini
date:       29/01/14
content:    Samples sequenced with PacBio
'''
# Modules
import pandas as pd



# Classes
class SampleTable(pd.DataFrame):
    def filter_seq_run(self, run):
        '''Get only the samples from one sequencing run'''
        return self[self.run == run]



# Globals
sample_list = [\
 {'name': 'S1', 'run': 'Upp23', 'description': 'NL4-3', 'filename': 'pb_023_1',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i')},
 {'name': 'S2', 'run': 'Upp23', 'description': 'Mix1', 'filename': 'pb_023_2',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i')},
 {'name': 'S3', 'run': 'Upp23', 'description': 'pat', 'filename': 'pb_023_3',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i')},
]
sample_table = SampleTable(sample_list)

# Make a dictionary
samples = {}
for line in sample_list:
    dic = dict(line)
    name = dic.pop('name')
    samples[name] = dic

