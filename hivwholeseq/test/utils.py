# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/01/15
content:    Support objects for tests
'''
# Classes
class Read(object):
    '''Mock read class'''
    qname = ''
    is_unmapped = False
    is_unpaired = True
    is_proper_pair = True
    is_reverse = False

    def __init__(self, seq, pos=0, **kwargs):
        self.seq = seq
        self.pos = pos

        if 'qual' not in kwargs:
            self.qual = 'G' * len(self.seq)

        if 'cigar' not in kwargs:
            self.cigar = [(0, len(self.seq))]

        for key, value in kwargs.iteritems():
            setattr(self, key, value)


    def __eq__(self, other):
        return (isinstance(other, self.__class__) and
                self.__dict__ == other.__dict__)


    def __str__(self):
        return self.__dict__.__str__()


    def __repr__(self):
        return self.__dict__.__repr__()



def fix_pair(reads):
    '''Fix insert size and mpos'''
    readf = reads[reads[0].is_reverse]
    readr = reads[reads[1].is_reverse]

    readf.mpos = readr.pos
    readr.mpos = readf.pos

    isize = (readr.pos - readf.pos +
             sum(bl for (bt, bl) in readr.cigar if bt in (0, 2)))

    readf.isize = isize
    readr.isize = -isize
