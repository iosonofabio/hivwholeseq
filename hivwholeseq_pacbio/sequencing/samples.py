# vim: fdm=marker
'''
author:     Fabio Zanini
date:       29/01/14
content:    Samples sequenced with PacBio
'''
# Modules
import numpy as np
import pandas as pd



# Classes
class SampleSeq(pd.Series):
    '''A sequenced sample (if something has been sequenced twice, they are separate)'''

    def __init__(self, *args, **kwargs):
        '''Initialize a sequenced sample'''
        super(SampleSeq, self).__init__(*args, **kwargs)
        self.name = self['name']

        from hivwholeseq_pacbio.sequencing.filenames import get_seqrun_foldername
        seq_run = self.loc['run']
        self['seqrun_folder'] = str(get_seqrun_foldername(seq_run))
        self['folder'] = self['seqrun_folder']+self['name']+'/'


    @property
    def _constructor(self):
        return SampleSeq


    def get_premapped_filename(self, **kwargs):
        '''Get the filename of the readed premapped to reference'''
        from hivwholeseq_pacbio.sequencing.filenames import get_premapped_filename as gfn
        return gfn(self.seqrun_folder, self.name, **kwargs)



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

 {'name': 'convPCR', 'run': 'Upp91', 'description': 'conventional PCR',
  'filename': 'pb_091_1', 'fragments': ['F4o']},
 {'name': 'emPCR', 'run': 'Upp91', 'description': 'emulsion PCR',
  'filename': 'pb_091_2', 'fragments': ['F4o']},
]
sample_table = SampleTable(sample_list)

# Make a dictionary
samples = {}
for line in sample_list:
    dic = dict(line)
    name = dic.pop('name')
    samples[name] = dic



# Functions
def load_sample(samplename):
    '''Load a sample'''
    sample  = samples[samplename]
    sample['name'] = samplename
    return SampleSeq(sample)
