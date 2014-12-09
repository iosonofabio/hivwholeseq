# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/12/14
content:    Support name for filenames.
'''
# Functions
def get_figure_folder(username, subfolder=''):
    if username == 'fzanini':
        fn = '/ebio/ag-neher/share/users/fzanini/phd/papers/hivintra/'
    elif username == 'fabio':
        fn = '/home/fabio/university/phd/papers/hivintra/'
    elif username in ('richard', 'rneher'):
        raise ValueError('R. please fill in this function in filenames.py')

    if subfolder:
        fn = fn + subfolder + '/'

    fn = fn + 'figures/'
    return fn
