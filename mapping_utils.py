# vim: fdm=indent
'''
author:     Fabio Zanini
date:       22/08/13
content:    Settings of stampy used by our mapping scripts.
'''
# Globals
stampy_bin = '/ebio/ag-neher/share/programs/bundles/stampy-1.0.22/stampy.py'
subsrate = '0.05'



# Functions
def get_ind_good_cigars(cigar, match_len_min=30):
    '''Keep only CIGAR blocks between two long matches'''
    from numpy import array

    criterion = lambda x: (x[0] == 0) and (x[1] >= match_len_min)
    good_cigars = array(map(criterion, cigar), bool, ndmin=1)

    # If there are no or one good CIGAR, keep that
    if (good_cigars).sum() < 2:
        return good_cigars
    
    # If there are 2+, keep stuff in the middle
    else:
        tmp = good_cigars.nonzero()[0]
        first_good_cigar = tmp[0]
        last_good_cigar = tmp[-1]
        good_cigars[first_good_cigar: last_good_cigar + 1] = True
    return good_cigars
