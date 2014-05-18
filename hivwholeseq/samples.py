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

 # First patient samples (nested PCR)
 {'name': 'VK04-3106_PCR2', 'run': 'Tue37', 'adaID': 'TS2', 'date': '2004-07-09',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i'),
  'patient': '20097'},
 {'name': '08HR-0235_PCR2', 'run': 'Tue37', 'adaID': 'TS4', 'date': '2008-02-07',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i'),
  'patient': '15823'},
 {'name': 'VK07-4778_PCR2', 'run': 'Tue37', 'adaID': 'TS5', 'date': '2007-07-12',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i'),
  'patient': '15313'},
 {'name': 'VK03-4298_PCR2', 'run': 'Tue37', 'adaID': 'TS6', 'date': '2003-09-30',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i'),
  'patient': '20097'},
 {'name': 'VK09-7738_PCR2', 'run': 'Tue37', 'adaID': 'TS7', 'date': '2009-09-24',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
  'patient': '20529'},
 {'name': 'VK08-8014_PCR2', 'run': 'Tue37', 'adaID': 'TS12', 'date': '2008-10-28',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
  'patient': '15376'},
 {'name': 'VK08-6634_PCR2', 'run': 'Tue37', 'adaID': 'TS13', 'date': '2008-09-11',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
  'patient': '20529'},
 {'name': 'VK07-8262_PCR2', 'run': 'Tue37', 'adaID': 'TS14', 'date': '2007-11-27',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
  'patient': '15376'},
 {'name': 'VK08-2987_PCR2', 'run': 'Tue37', 'adaID': 'TS15', 'date': '2008-04-21',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
  'patient': '15313'},

 # Test samples for PCR1, Nextera XT, and MiSeq v3 (600 bp)
 {'name': 'MIX1_PCR1', 'run': 'Tue42', 'adaID': 'N1-S1',
  'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5bo', 'F6o'),
  'description': 'MIX1 SF162 50% + NL4-3 50%'},
 {'name': 'MIX2_PCR1', 'run': 'Tue42', 'adaID': 'N3-S3',
  'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5bo', 'F6o'),
  'description': 'MIX2 SF162 95% + NL4-3 4.5% + F10 0.5%'},
 {'name': 'NL4-3_PCR1_TaqHiFi', 'run': 'Tue42', 'adaID': 'N4-S3',
  'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5bo', 'F6o')},
 {'name': 'NL4-3_PCR1_Taq', 'run': 'Tue42', 'adaID': 'N5-S4',
  'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5bo', 'F6o')},
 {'name': 'NL4-3_PCR2_TaqHiFi', 'run': 'Tue42', 'adaID': 'N2-S2',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i')},
 {'name': 'NL4-3_PCR2_Taq', 'run': 'Tue42', 'adaID': 'N6-S4',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5bi', 'F6i')},

 # TrueSeq Nano patient samples #2
 {'name': 'F10_PCR2_TaqHiFi', 'run': 'Tue44', 'adaID': 'TS2',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i')},
 {'name': 'F10_PCR2_Taq', 'run': 'Tue44', 'adaID': 'TS4',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i')},
 {'name': '34493_PCR2', 'run': 'Tue44', 'adaID': 'TS5',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
  'patient': '15107', 'date': '2006-05-16'},
 {'name': '29184_PCR2', 'run': 'Tue44', 'adaID': 'TS6',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
  'patient': '9669', 'date': '2004-05-05'},
 {'name': '30847_PCR2', 'run': 'Tue44', 'adaID': 'TS7',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
  'patient': '15241', 'date': '2004-12-27'},
 {'name': 'VK02-4452_PCR2', 'run': 'Tue44', 'adaID': 'TS12',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
  'patient': '20097', 'date': '2002-10-29'},
 {'name': '06HR-0145_PCR2', 'run': 'Tue44', 'adaID': 'TS13',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
  'patient': '15823', 'date': '2006-02-02'},
 {'name': '28929_PCR2', 'run': 'Tue44', 'adaID': 'TS14',
  'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
  'patient': '15107', 'date': '2004-03-30'},

 # Nextera XT #2 + BluePippin test + 3 pat PCR1 samples
 {'name': 'F10_PCR1_TaqHiFi', 'run': 'Tue48', 'adaID': 'N2-S2',
  'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o')},
 {'name': 'F10_PCR1_Taq', 'run': 'Tue48', 'adaID': 'N3-S2',
  'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o')},
 {'name': 'RNA_mix_PCR1_TaqHiFi', 'run': 'Tue48', 'adaID': 'N4-S1',
  'fragments': ('F2o', 'F3o', 'F4o', 'F6o')},
 {'name': 'RNA_mix_PCR1_Taq', 'run': 'Tue48', 'adaID': 'N5-S1',
  'fragments': ('F2o', 'F3o', 'F4o', 'F6o')},
 {'name': 'RNA_mix_PCR2_TaqHiFi', 'run': 'Tue48', 'adaID': 'N6-S1',
  'fragments': ('F2i', 'F3i', 'F4i', 'F6i')},
 {'name': 'RNA_mix_PCR2_Taq', 'run': 'Tue48', 'adaID': 'N1-S3',
  'fragments': ('F2i', 'F3i', 'F4i', 'F6i')},
 {'name': 'VK02-4452_PCR1', 'run': 'Tue48', 'adaID': 'N2-S3',
  'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o'),
  'patient': '20097', 'date': '2002-10-29'},
 {'name': '06HR-0145_PCR1', 'run': 'Tue48', 'adaID': 'N1-S4',
  'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o'),
  'patient': '15823', 'date': '2006-02-02'},
 {'name': '28929_PCR1', 'run': 'Tue48', 'adaID': 'N2-S4',
  'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o')},

  ## TruSeq #3 + South Africa #2 (mixed library, only read1 worked)
  #{'name': 'POL_SA_2', 'run': 'Tue52', 'adaID': 'TS1',
  # 'description': 'South Africa mutation rate #2'}, #SOUTH AFRICAN SAMPLE
  #{'name': 'VK01-2965_PCR2', 'run': 'Tue52', 'adaID': 'TS5',
  # 'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'), #IS THIS TRUE?
  # 'patient': '20097', 'date': '2001-08-13'},
  #{'name': '05HR-0269_PCR2', 'run': 'Tue52', 'adaID': 'TS6',
  # 'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'), #IS THIS TRUE?
  # 'patient': '15823', 'date': '2005-03-01'},
  #{'name': 'VK09-1685_PCR2', 'run': 'Tue52', 'adaID': 'TS15',
  # 'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'), #IS THIS TRUE?
  # 'patient': '15313', 'date': '2009-02-26'},

  # reference RNA and patients PCR1 reruns
  {'name': 'LAI-III', 'run': 'Tue59', 'adaID': 'N3-S3',
   'description': 'Reference sequence for subtype B',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o')},
  {'name': '38540', 'run': 'Tue59', 'adaID': 'N5-S4',
   'description': 'Reference sequence for subtype C',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o')},
  {'name': '38304', 'run': 'Tue59', 'adaID': 'N6-S3',
   'description': 'Reference sequence for subtype B',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o')},
  {'name': 'VK03-4298_PCR1', 'run': 'Tue59', 'adaID': 'N5-S2',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5bo', 'F6o'),
   'patient': '20097',
   'date': '2003-09-30'},
  {'name': 'VK04-3106_PCR1', 'run': 'Tue59', 'adaID': 'N3-S4',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5bo', 'F6o'),
   'patient': '20097',
   'date': '2004-07-09'},
  {'name': 'VK09-1685_PCR1', 'run': 'Tue59', 'adaID': 'N4-S2',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o'),
   'patient': '15313', 'date': '2009-02-26'},
  {'name': 'VK09-7738_PCR1', 'run': 'Tue59', 'adaID': 'N6-S4',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o'),
   'patient': '20529', 'date': '2009-09-24'},

  # 14 patient samples with Nextera XT
  {'name': '28338_PCR1', 'run': 'Tuen3', 'adaID': 'N1-S1',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o'),
   'patient': '15241', 'date': '2004-01-15'},
  {'name': '29698_PCR1', 'run': 'Tuen3', 'adaID': 'N2-S1',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o'),
   'patient': '15319', 'date': '2004-07-19'},
  {'name': '21006_PCR1', 'run': 'Tuen3', 'adaID': 'N3-S1',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o'),
   'patient': '9669', 'date': '2001-01-16'},
  {'name': 'VK08-1001_PCR1', 'run': 'Tuen3', 'adaID': 'N4-S1',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o'),
   'patient': '15363', 'date': '2008-02-05'},
  {'name': '24890_PCR1', 'run': 'Tuen3', 'adaID': 'N5-S2',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o'),
   'patient': '15107', 'date': '2002-08-19'},
  {'name': '04HR-1501_PCR1', 'run': 'Tuen3', 'adaID': 'N6-S2',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o'),
   'patient': '15823', 'date': '2004-12-19'},
  {'name': 'VK06-1885_PCR1', 'run': 'Tuen3', 'adaID': 'N1-S2',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o'),
   'patient': '20529', 'date': '2006-04-10'},
  {'name': '25775_PCR1', 'run': 'Tuen3', 'adaID': 'N2-S2',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o'),
   'patient': '15241', 'date': '2002-12-17'},
  {'name': 'VK06-6001_PCR1', 'run': 'Tuen3', 'adaID': 'N3-S3',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o'),
   'patient': '15376', 'date': '2006-10-24'},
  {'name': '21484_PCR1', 'run': 'Tuen3', 'adaID': 'N4-S3',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o'),
   'patient': '15107', 'date': '2001-03-26'},
  {'name': '26477_PCR1', 'run': 'Tuen3', 'adaID': 'N5-S3',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o'),
   'patient': '9669', 'date': '2003-04-03'},
  {'name': 'VK00-1524_PCR1', 'run': 'Tuen3', 'adaID': 'N6-S4',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o'),
   'patient': '20097', 'date': '2000-05-02'},
  {'name': 'VK07-4218_PCR1', 'run': 'Tuen3', 'adaID': 'N1-S4',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o'),
   'patient': '20529', 'date': '2007-06-15'},
  {'name': '7686_PCR1', 'run': 'Tuen3', 'adaID': 'N2-S4',
   'fragments': ('F1o', 'F2o', 'F3o', 'F4o', 'F5ao', 'F6o'),
   'patient': '15034', 'date': '1994-08-25'},

  # emPCR, SA#3 plasmid, and SE#5 PCR2
  {'name': 'POL_SA_3_plasmid', 'run': 'Tuen6', 'adaID': 'TS1',
   'description': 'South Africa mutation rate #3, plasmid control'}, #SOUTH AFRICAN SAMPLE
  {'name': 'convPCR1542', 'run': 'Tuen6', 'adaID': 'TS2',
   'description': 'emulsion PCR test: conventional PCR control with 1542 bp long amplicon'},
  {'name': 'convPCR1968', 'run': 'Tuen6', 'adaID': 'TS5',
   'description': 'emulsion PCR test: conventional PCR control with 1968 bp long amplicon'},
  {'name': 'emPCR1542', 'run': 'Tuen6', 'adaID': 'TS6',
   'description': 'emulsion PCR test: emulsion PCR with 1542 bp long amplicon'},
  {'name': 'emPCR1968', 'run': 'Tuen6', 'adaID': 'TS12',
   'description': 'emulsion PCR test: emulsion PCR with 1968 bp long amplicon'},
  {'name': '28338_PCR2', 'run': 'Tuen6', 'adaID': 'TS18',
   'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
   'patient': '15241', 'date': '2004-01-15'},
  {'name': '29698_PCR2', 'run': 'Tuen6', 'adaID': 'TS19',
   'fragments': ('F1i', 'F2i', 'F3i', 'F4i', 'F5ai', 'F6i'),
   'patient': '15319', 'date': '2004-07-19'},

 ]

# Sample table using Pandas
sample_table = SampleTable(sample_list)

# Make a dictionary
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
