# vim: fdm=indent
'''
author:     Fabio Zanini
date:       12/03/15
content:    Support functions for pandas.
'''
# Modules
from __future__ import absolute_import
import numpy as np
import pandas as pd



# Functions
def add_binned_column(data, name, column, bins=10, clip=False):
    '''Bin the contents of a column and add it as an additional column'''
    if np.isscalar(bins):
        bins = np.unique(np.array(data.loc[:, column].quantile(q=np.linspace(0, 1, bins + 1))))
    data[name] = pd.cut(data.loc[:, column], bins=bins, include_lowest=True, labels=False)
    if clip:
        data[name] = data[name].clip(0, len(bins) - 2)

    binsc = 0.5 * (bins[1:] + bins[:-1])
    return bins, binsc
