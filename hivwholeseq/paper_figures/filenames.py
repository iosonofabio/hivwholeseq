# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/12/14
content:    Support name for filenames.
'''
# Functions
def get_root_folder(username, subfolder=''):
    if username == 'fzanini':
        fn = '/ebio/ag-neher/share/users/fzanini/phd/papers/'
    elif username == 'fabio':
        fn = '/home/fabio/university/phd/papers/'
    elif username in ('richard', 'rneher'):
        fn = '/ebio/ag-neher/share/users/rneher/'

    if subfolder in ['first', 'controls']:
        fn = fn + 'HIVEVO_first_figures/'

    elif subfolder == 'popgen':
        fn = fn + 'HIVEVO_popgen/'

    return fn


def get_figure_folder(username, subfolder=''):
    return get_root_folder(username, subfolder=subfolder)+'figures/'


def get_data_folder(username, subfolder=''):
    return get_root_folder(username, subfolder=subfolder)+'data/'
