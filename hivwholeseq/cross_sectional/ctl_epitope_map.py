# vim: fdm=indent
'''
author:     Fabio Zanini
date:       03/05/15
content:    Parse the epitope map from LANL.
'''
# Modules



# Functions
def get_ctl_epitope_map(species=None):
    '''Get the CTL epitope map from LANL'''
    import pandas as pd
    from hivwholeseq.cross_sectional.filenames import (
        get_ctl_epitope_map_filename)
    fn = get_ctl_epitope_map_filename()
    table = pd.read_csv(fn, skiprows=1)

    if species is not None:
        table = table.loc[table['Species'] == species]
        del table['Species']

    table['Protein'] = [x.lower() for x in table['Protein']]

    return table



# Script
if __name__ == '__main__':

    table = get_ctl_epitope_map()
