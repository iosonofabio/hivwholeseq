# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/12/14
content:    Support name for filenames.
'''
# Functions
def get_figure_folder(username, subfolder=''):
    if username == 'fzanini':
        fn = '/ebio/ag-neher/share/users/fzanini/phd/papers/'
    elif username == 'fabio':
        fn = '/home/fabio/university/phd/papers/'
    elif username in ('richard', 'rneher'):
        fn = '/ebio/ag-neher/share/users/rneher/'

    if subfolder in ['controls']:
        fn = fn + 'HIVEVO_support/' +subfolder + '/figures/'

    elif subfolder == 'first':
        fn = fn + 'HIVEVO_first_paper/figures/'

    elif subfolder == 'popgen':
        fn = fn + 'HIVEVO_popgen/figures/'

    return fn
