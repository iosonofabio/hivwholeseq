# vim: fdm=marker
'''
author:     Fabio Zanini
date:       29/01/14
content:    Samples sequenced with PacBio
'''
# Modules
import pandas as pd


# Globals
samples = [{'name': 'S1', 'run': 'Upp23', 'description': 'NL4-3', 'filename': 'pb_023_1'},
           {'name': 'S2', 'run': 'Upp23', 'description': 'Mix1', 'filename': 'pb_023_2'},
           {'name': 'S3', 'run': 'Upp23', 'description': 'pat', 'filename': 'pb_023_3'},
           ]
samples = pd.DataFrame(samples)

