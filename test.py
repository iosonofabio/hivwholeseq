def pair_generator(iterable):
    it = iter(iterable)
    while True:
        try:
            a = it.next()
            b = it.next()
            yield (a, b)
        except StopIteration:
            raise
#
#a = [1, 2, 3, 4]
#b = [1, 2, 3, 4, 5]
#
#print [aa for aa in pair_generator(a)]
#print [aa for aa in pair_generator(b)]

def identity(x):
    '''Test function for import only'''
    return x


print __file__


bs = []
for aa in a:
    b = np.zeros((4, 6, len(aa[0][0])), int)
    for i, aaa in enumerate(aa):
        for ii, aaaa in enumerate(aaa):
            b[i, ii] = aaaa
    bs.append(b)
