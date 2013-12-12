# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/08/13
content:    Description module for our HIV datasets.
'''
# Globals
test_tiny = {'description': 'Small subsample of run 37 for testing purposes',
             'date': '2013-12-06',
             'n_cycles': 500,
             'comments': 'MiSeq run37 (tiny subsample)',
             'folder': '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/test_tiny/',
             'raw_data': {'read1': 'read1.fastq',
                          'adapter': 'barcode.fastq',
                          'read2': 'read2.fastq'},
             'adapters': ('TS2',  'TS4',  'TS5',  'TS6', 'TS7',
                          'TS12', 'TS13', 'TS14', 'TS15'),
             'samples': ('VK04-3106-test',
                         '08HR-0235-test',
                         'VK07-4778-test',
                         'VK03-4298-test',
                         'VK09-7738-test',
                         'VK08-8014-test',
                         'VK08-6634-test',
                         'VK07-8262-test',
                         'VK08-2987-test',
                        ),
            }


testmiseq = {'description': 'Test run for the MiSeq and the PCR',
             'date': '2013-07-30',
             'n_cycles': 500,
             'comments': 'MiSeq run28',
             'folder': '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/',
             'raw_data': {'read1': 'lane1_NoIndex_L001_R1_001.fastq',
                          'adapter': 'lane1_NoIndex_L001_R2_001.fastq',
                          'read2': 'lane1_NoIndex_L001_R3_001.fastq'},
             'adapters': ('TS2', 'TS4', 'TS7', 'TS16', 'TS18', 'TS19'),
             'samples': ('NL4-3',
                         'SF162',
                         'F10',
                         '37024',
                         'MIX1',
                         'MIX2',
                        ),
            }

patients1 = {'description': 'Second library (first patient-only)',
             'date': '2013-09-30',
             'n_cycles': 500,
             'comments': 'MiSeq run37',
             'folder': '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run37/',
             'raw_data': {'read1': 'read1.fastq',
                          'adapter': 'barcode.fastq',
                          'read2': 'read2.fastq'},
             'adapters': ('TS2', 'TS4', 'TS5', 'TS6', 'TS7',
                          'TS12', 'TS13', 'TS14', 'TS15'),
             'samples': ('VK04-3106',
                         '08HR-0235',
                         'VK07-4778',
                         'VK03-4298',
                         'VK09-7738',
                         'VK08-8014',
                         'VK08-6634',
                         'VK07-8262',
                         'VK08-2987',
                        ),
            }


testnextera_Lina = {'description': 'Test run for the Nextera library (Sweden)',
                    'date': '2013-10-20',
                    'n_cycles': 300,
                    'comments': 'Library prepared by Lina, sequencing by Pawel Zajac (illumina)',
                    'folder': '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/pawel_nextera/',
                    'raw_data': None, # ALREADY DEMULTIPLEXED!
                    'adapters': ('01', '02'), # FAKE NUMBERS!
                    'samples': ('Nextera_HIV-8262-1',
                                'Nextera_HIV-8262-2'),
                   }


# MiSeq runs
MiSeq_runs = {'Tue28': testmiseq,
              'Tue37': patients1,
              'Lina_nextera': testnextera_Lina,
              'test_tiny': test_tiny,
             }
