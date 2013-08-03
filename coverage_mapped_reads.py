# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/08/13
content:    Check the coverage of mapped reads on the HIV genome.
'''
# Modules
import os
import sys
from Bio import SeqIO
import numpy as np
import pysam
from map_HIV_HXB2 import load_adapter_table



# Globals
VERBOSE = 1
data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/copies/'
adapters_table_file = 'adapters_table.dat'
mapped_file = 'mapped_to_HXB2.sam'



# Script
if __name__ == '__main__':

    # Directory to read
    adapter_table = load_adapter_table(data_folder)
    adaID = adapter_table['ID'][adapter_table['sample'] == 'NL4-3'][0]
    dirname = 'adapterID_'+'{:02d}'.format(adaID)+'/'

    # SAM file (output of stampy)
    samfilename = data_folder+dirname+mapped_file
    samfile = pysam.Samfile(samfilename, 'r')


