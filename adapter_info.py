# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/08/13
content:    Info on the illumina adapters.
'''
# Globals
# From: http://support.illumina.com/downloads/illumina_adapter_sequences_letter.ilmn
adapters_LT = {1: 'ATCACG',
               2: 'CGATGT',
               3: 'TTAGGC',
               4: 'TGACCA',
               5: 'ACAGTG',
               6: 'GCCAAT',
               7: 'CAGATC',
               8: 'ACTTGA',
               9: 'GATCAG',
               10: 'TAGCTT',
               11: 'GGCTAC',
               12: 'CTTGTA',
               13: 'AGTCAA',
               14: 'AGTTCC',
               15: 'ATGTCA',
               16: 'CCGTCC',
               18: 'GTCCGC',
               19: 'GTGAAA',
               20: 'GTGGCC',
               21: 'GTTTCG',
               22: 'CGTACG',
               23: 'GAGTGG',
               25: 'ACTGAT',
               27: 'ATTCCT'}


adapters_table_file = 'adapters_table.dat'




# Functions
def load_adapter_table(data_folder):
    '''Load table of adapters and samples'''
    from numpy import loadtxt
    table = loadtxt(data_folder+adapters_table_file,
                    dtype=[('seq', 'S6'), ('ID', int), ('sample', 'S50')],
                    ndmin=1)
    return table


def foldername_adapter(adaID):
    '''Convert an adapter number in a folder name'''
    return 'adapterID_'+'{:02d}'.format(adaID)+'/'
