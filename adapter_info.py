# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/08/13
content:    Info on the illumina adapters.
'''
# Globals
# From: http://support.illumina.com/downloads/illumina_adapter_sequences_letter.ilmn
adapters_LT = {01: 'ATCACG',
               02: 'CGATGT',
               03: 'TTAGGC',
               04: 'TGACCA',
               05: 'ACAGTG',
               06: 'GCCAAT',
               07: 'CAGATC',
               08: 'ACTTGA',
               09: 'GATCAG',
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
def foldername_adapter(adaID):
    '''Convert an adapter number in a folder name'''
    return 'adapterID_'+'{:02d}'.format(adaID)+'/'
