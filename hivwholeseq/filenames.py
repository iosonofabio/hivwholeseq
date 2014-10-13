# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/10/14
content:    Module containing basic paths for the whole analysis.
'''
# Modules
import os

import hivwholeseq


# Globals
# FIXME: use env vars and similia
if os.path.isdir('/media/FZ_MPI/HIV_Sweden/'):
    root_data_folder = '/media/FZ_MPI/HIV_Sweden/'
    stampy_bin = '/usr/bin/stampy'
    spades_bin = None
    bwa_bin = None
else:
    root_data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/'
    stampy_bin = '/ebio/ag-neher/share/programs/bundles/stampy-1.0.22/stampy.py'
    bwa_bin = '/ebio/ag-neher/share/programs/bin/bwa'
    spades_bin = '/ebio/ag-neher/share/programs/bundles/SPAdes-2.5.0-Linux/bin/spades.py'


reference_folder = root_data_folder+'reference/'
table_filename = hivwholeseq.__path__[0] + '/tables/HIV_reservoir_all.xlsx'
