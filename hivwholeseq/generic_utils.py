# vim: fdm=marker
'''
author:     Fabio Zanini
date:       25/10/13
content:    A bunch of generic utility functions that should be in Python but are
            not.
'''
# Functions 
def mkdirs(newdir):
    """ Create a directory and all parent folders.
        Features:
        - parent directoryies will be created
        - if directory already exists, then do nothing
        - if there is another filsystem object with the same name, raise an exception
    """
    import os
    if os.path.isdir(newdir):
        pass
    elif os.path.isfile(newdir):
        raise OSError("cannot create directory, file already exists: '%s'" % newdir)
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            mkdirs(head)
        if tail:
            os.mkdir(newdir)

