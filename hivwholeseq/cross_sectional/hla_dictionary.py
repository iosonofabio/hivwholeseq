# vim: fdm=indent
'''
author:     Fabio Zanini
date:       04/05/15
content:    Load the HLA dictionary to map between genotypes, subtypes,
            serotypes, supertypes, and other madnesses of that kind.
'''
# Modules



# Functions
def get_hla_dictionary(locus=None):
    '''Get the dictionary of HLA to serotype/supertype'''
    import pandas as pd
    from hivwholeseq.cross_sectional.filenames import (
        get_hla_dictionary_filename)
    fn = get_hla_dictionary_filename()
    table = pd.read_csv(fn, skiprows=0)

    if locus is not None:
        if locus in ('A', 'B', 'C'):
            table = table.loc[[x[0] == locus for x in table['Genotype']]]

    return table



# Script
if __name__ == '__main__':

    table = get_hla_dictionary()
