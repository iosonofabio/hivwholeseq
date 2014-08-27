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


def getchar():
    '''Get single char from terminal, much like getchar() in C'''
    import os
    import sys
    import termios

    fd = sys.stdin.fileno()
    if os.isatty(fd):
    	old = termios.tcgetattr(fd)
    	new = termios.tcgetattr(fd)
    	new[3] = new[3] & ~termios.ICANON & ~termios.ECHO
    	new[6] [termios.VMIN] = 1
    	new[6] [termios.VTIME] = 0

    	try:
    	    termios.tcsetattr(fd, termios.TCSANOW, new)
    	    termios.tcsendbreak(fd,0)
    	    ch = os.read(fd,7)
    	finally:
    	    termios.tcsetattr(fd, termios.TCSAFLUSH, old)
    else:
        ch = os.read(fd,7)
    
    return(ch)
