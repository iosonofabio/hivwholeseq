# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/08/13
content:    Description module for our HIV datasets.
'''
# Modules
from hivwholeseq.samples import sample_table


# Globals
MiSeq_runs_list = [\
 {'name': 'test_tiny',
  'description': 'Small subsample of run 37 for testing purposes',
  'date': '2013-12-06',
  'n_cycles': 500,
  'comments': 'MiSeq run37 (tiny subsample)',
  'folder': '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/test_tiny/',
  'raw_data': {'read1': 'read1.fastq',
               'adapter': 'barcode.fastq',
               'read2': 'read2.fastq'},
 },

 {'name': 'Tue28',
  'description': 'Test run for the MiSeq and the PCR',
  'date': '2013-07-30',
  'n_cycles': 500,
  'comments': 'MiSeq run28',
  'folder': '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/',
  'raw_data': {'read1': 'lane1_NoIndex_L001_R1_001.fastq',
               'adapter': 'lane1_NoIndex_L001_R2_001.fastq',
               'read2': 'lane1_NoIndex_L001_R3_001.fastq'},
 },

 {'name': 'Tue37',
  'description': 'Second library (first patient-only)',
  'date': '2013-09-30',
  'n_cycles': 500,
  'comments': 'MiSeq run37',
  'folder': '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run37/',
  'raw_data': {'read1': 'read1.fastq',
               'adapter': 'barcode.fastq',
               'read2': 'read2.fastq'},
 },

 {'name': 'Lina_nextera',
  'description': 'Test run for the Nextera library (Sweden)',
  'date': '2013-10-20',
  'n_cycles': 300,
  'comments': 'Library prepared by Lina, sequencing by Pawel Zajac (illumina)',
  'folder': '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/pawel_nextera/',
  'raw_data': None, # ALREADY DEMULTIPLEXED!
 },

 {'name': 'Tue42',
  'description': 'Test run for PCR1 + Nextera XT + 600 bp',
  'date': '2014-01-10',
  'n_cycles': 602,
  'comments': 'Dual index Nextera adapters',
  'folder': '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run42/',
  'raw_data': {'read1': 'lane1_NoIndex_L001_R1_001.fastq',
               'read2': 'lane1_NoIndex_L001_R4_001.fastq',
               'adapter1': 'lane1_NoIndex_L001_R2_001.fastq',
               'adapter2': 'lane1_NoIndex_L001_R3_001.fastq'},
 },

 {'name': 'Tue44',
  'description': 'Second patient-only library with TrueSeq Nano',
  'date': '2014-01-24',
  'n_cycles': 602,
  'folder': '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run44/',
  'raw_data': {'read1': 'lane1_NoIndex_L001_R1_001.fastq',
               'read2': 'lane1_NoIndex_L001_R3_001.fastq',
               'adapter': 'lane1_NoIndex_L001_R2_001.fastq'},
 },

 {'name': 'Tue48',
  'description': 'Second Nextera XT + test for BluePippin',
  'date': '2014-01-31',
  'n_cycles': 600,
  'folder': '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run48/',
  'raw_data': None, # ALREADY DEMULTIPLEXED!
 },

 {'name': 'Tue52',
  'description': 'Third library with TrueSeq Nano',
  'date': '2014-02-18',
  'n_cycles': 600,
  'folder': '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run52/',
  'raw_data': {'read1': 'lane1_NoIndex_L001_R1_001.fastq.gz',
               'read2': 'lane1_NoIndex_L001_R3_001.fastq.gz',
               'adapter': 'lane1_NoIndex_L001_R2_001.fastq.gz'},
 },


]

# Add samples to the runs
for run in MiSeq_runs_list:
    tmp = sample_table.filter_seq_run(run['name'])
    run['adapters'] = tuple(tmp['adaID'])
    run['samples'] = tuple(tmp['name'])
    run['sample_table'] = tmp.copy()

# runs dictionary
MiSeq_runs = {run['name']: run for run in MiSeq_runs_list}
