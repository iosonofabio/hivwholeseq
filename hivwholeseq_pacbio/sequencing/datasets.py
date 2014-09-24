# vim: fdm=marker
'''
author:     Fabio Zanini
date:       30/01/14
content:    Datasets (sequencing runs)
'''
# Modules
from hivwholeseq_pacbio.sequencing.samples import sample_table
from hivwholeseq_pacbio.sequencing.filenames import get_seqrun_foldername



# Globals
PacBio_runs_list = [\
 
 {'name': 'Upp23',
  'description': 'First test run',
  'date': '2014-01-21'},

 {'name': 'Upp91',
  'description': 'Second test run (emPCR)',
  'date': '2014-09-22'},

]

# Add samples to the runs
for run in PacBio_runs_list:
    tmp = sample_table.filter_seq_run(run['name'])
    run['samples'] = tuple(tmp['name'])
    run['sample_table'] = tmp.copy()
    run['folder'] = get_seqrun_foldername(run['name'])

# runs dictionary
PacBio_runs = {run['name']: run for run in PacBio_runs_list}
