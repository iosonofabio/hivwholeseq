# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/01/15
content:    Support objects for tests
'''
# Classes
class Read(object):
    '''Mock read class'''
    def __eq__(self, other):
        return (isinstance(other, self.__class__) and
                self.__dict__ == other.__dict__)

