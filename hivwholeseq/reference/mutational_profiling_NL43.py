# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/01/15
content:    Al-Mawsawi et al (Retrovirology 2014) have profiled for fitness many
            single-mutants of NL4-3. I collect those data here.
'''
# Modules
import os
import sys
import numpy as np

from hivwholeseq.reference import load_custom_reference


# Functions
def load_table(align_to_reference=None):
    '''Load the TSV table for the mutational profile'''
    import pandas as pd

    from hivwholeseq.filenames import table_folder
    table_filename = table_folder+'Al-Mawsawi_etal_2014/s12977-014-0124-6-s6-sheet1.tsv'

    table = pd.read_csv(table_filename,
                        header=6,
                        sep='\t',
                        usecols=(0, 1, 5),
                        dtype={'Protein': 'S10',
                               'Mutation': 'S10',
                               'RC Index': float},
                       )

    table['Pos NL4-3'] = table['Mutation'].str[1:-1].astype(int) - 1
    table['Ancestral'] = table['Mutation'].str[0]
    table['Derived'] = table['Mutation'].str[-1]

    if align_to_reference is not None:
        table = add_coordinates_reference(table, refname=align_to_reference)

    return table
    

def add_coordinates_reference(table, refname='HXB2', filter_missing=True, VERBOSE=0):
    '''Add coordinates of another reference'''
    from seqanpy import align_global

    ref1 = load_custom_reference('NL4-3')
    ref2 = load_custom_reference(refname)

    (score, ali1, ali2) = align_global(ref1, ref2, score_gapopen=-20)

    alim1 = np.fromstring(ali1, 'S1')
    alim2 = np.fromstring(ali2, 'S1')
    poss = np.vstack([(alim1 != '-').cumsum() - 1,
                      (alim2 != '-').cumsum() - 1])

    # Exclude gaps
    if filter_missing:
        ind = (alim1 != '-') & (alim2 != '-')
    else:
        ind = (alim1 != '-')
    poss = poss[:, ind]

    # Make dictionary for quick access
    posdict = dict(poss.T)

    # Fill the new vector
    possref1 = table['Pos NL4-3']
    possref2 = - np.ones(len(possref1), int)
    for i, pos1 in enumerate(possref1):
        if pos1 not in posdict:
            continue

        pos2 = posdict[pos1]
        possref2[i] = pos2

        if VERBOSE >= 3:
            print ref1[pos1], ref2[pos2]

    # Add to table
    table['Pos '+refname] = possref2

    if filter_missing:
        table = table.loc[table['Pos '+refname] != -1]

    return table



# Script
if __name__ == '__main__':

    table = load_table()

    # Check that the ancestral alleles correspond to NL4-3, i.e. figure out what
    # the hell those coordinates are (NL4-3? HXB2?)
    print 'Checking coordinates with NL4-3...',
    ref = load_custom_reference('NL4-3')
    for _, row in table.iterrows():
        pos = row['Pos NL4-3']
        anc = row['Ancestral']
        if anc != ref[pos]:
            print pos, anc, ref[pos]
            break
    else:
        print 'OK'


    table = add_coordinates_reference(table, VERBOSE=2)
