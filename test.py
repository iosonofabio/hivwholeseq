def pair_generator(iterable):
    it = iter(iterable)
    while True:
        try:
            a = it.next()
            b = it.next()
            yield (a, b)
        except StopIteration:
            raise

a = [1, 2, 3, 4]
b = [1, 2, 3, 4, 5]

print [aa for aa in pair_generator(a)]
print [aa for aa in pair_generator(b)]
