# vim: fdm=indent
'''
author:     Fabio Zanini
date:       03/05/15
content:    Parse the epitope map from LANL.
'''
# Modules



# Functions
def get_ctl_epitope_map(species='human', hlas=None):
    '''Get the CTL epitope map from LANL
    
    Parameters:
       species (str): limit to epitopes for a certain species (e.g. human)
       hlas (list): limit to epitopes specific to any of these HLA
    '''
    import pandas as pd
    from hivwholeseq.cross_sectional.filenames import (
        get_ctl_epitope_map_filename)
    fn = get_ctl_epitope_map_filename()
    table = pd.read_csv(fn, skiprows=1)

    if species is not None:
        table = table.loc[table['Species'] == species]
        del table['Species']

    if hlas is not None:
       pass 

    table['Protein'] = [x.lower() for x in table['Protein']]

    return table


def extend_hla(hla):
    '''Extend HLA of patients to generic ones'''
    hla_ext = []
    for item in hla:
        item = item.replace(':', '')
        hla_ext.append(item)
        hla_ext.append(item[:6])
        hla_ext.append(item[:4])
    
    return hla_ext


def get_ctl_epitope_hla(table, hla):
    '''Get CTL epitopes specific to some HLA (patient)'''
    hla_ext = extend_hla(hla)
    ind = [i for i, item in enumerate(table['HLA'])
           if any(h in str(item) for h in hla_ext)]

    return table.iloc[ind]




# Script
if __name__ == '__main__':

    table = get_ctl_epitope_map()
