# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/08/13
content:    Description module for our HIV datasets.
'''
# Globals
dataset_testmiseq = {'description': 'Test run for the MiSeq and the PCR',
                     'adapters': (2, 4, 7, 16, 18, 19),
                     'samples': ('NL4-3', 'SF162', 'F10',
                                 'patient 37024',
                                 'MIX1 SF162 50% + NL4-3 50%',
                                 'MIX2 SF162 95% + NL4-3 4.5% + F10 0.5%'),
                     'folder': '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/',
                     'date': '2013-07-30',
                     'comments': 'MiSeq run28',
                     'primerF5': tuple('F5b' for s in xrange(6)),
                     'raw_data': {'read1': 'lane1_NoIndex_L001_R1_001.fastq',
                                  'adapter': 'lane1_NoIndex_L001_R2_001.fastq',
                                  'read2': 'lane1_NoIndex_L001_R3_001.fastq'},
                     'n_cycles': 500,
                    }

dataset_2 = {'description': 'Second library (first patient-only)',
             'adapters': (2, 4, 5, 6, 7, 12, 13, 14, 15),
             'samples': ('VK04-3106', '08HR-0235', 'VK07-4778', 'VK03-4298',
                         'VK09-7738', 'VK08-8014', 'VK08-6634', 'VK07-8262',
                         'VK08-2987'),
             'folder': '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run37/',
             'date': '2013-09-30',
             'comments': 'MiSeq run37',
             'primerF5': ('F5b', 'F5b', 'F5b', 'F5b',
                          'F5a', 'F5a', 'F5a', 'F5a', 'F5a'),
             'raw_data': {'read1': 'read1.fastq',
                          'adapter': 'barcode.fastq',
                          'read2': 'read2.fastq'},
             'n_cycles': 500,
            }


dataset_nextera = {'description': 'Test run for the Nextera library (Sweden)',
                   'adapters': (1, 2),
                   'samples': ('HIV-8262-1: 0.2 ng/ul', 'HIV-8262-2: 0.033 ng/ul'),
                   'folder': '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/pawel_nextera/',
                   'date': '2013-10-20',
                   'comments': 'Library prepared by Lina, sequencing by Pawel Zajac (illumina)',
                   'primerF5': ('F5a', 'F5a'), # FIXME: Find out about this
                   'n_cycles': 300,
                  }


# MiSeq runs
MiSeq_runs = {28: dataset_testmiseq,
              37: dataset_2,
              -1: dataset_nextera,
             }
