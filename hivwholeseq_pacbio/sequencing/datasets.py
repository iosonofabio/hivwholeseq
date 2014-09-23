# vim: fdm=marker
'''
author:     Fabio Zanini
date:       30/01/14
content:    Datasets (sequencing runs)
'''
# Modules
from hivwholeseq_pacbio.sequencing.samples import sample_table



# Globals
PacBio_runs_list = [\
 
 {'name': 'Upp23',
  'description': 'First test run',
  'date': '2014-01-21',
  'folder': '/ebio/ag-neher/share/data/PacBio_HIV_Karolinska/run23/'},

 {'name': 'Upp91',
  'description': 'Second test run (emPCR)',
  'date': '2014-09-22',
  'folder': '/ebio/ag-neher/share/data/PacBio_HIV_Karolinska/run91/'},

]

# Add samples to the runs
for run in PacBio_runs_list:
    tmp = sample_table.filter_seq_run(run['name'])
    run['samples'] = tuple(tmp['name'])
    run['sample_table'] = tmp.copy()

# runs dictionary
PacBio_runs = {run['name']: run for run in PacBio_runs_list}
