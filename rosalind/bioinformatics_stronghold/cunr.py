# Counting Unrooted Binary Trees

from functools import reduce


# The answer here is just (2n -5)!!
# We can help the calculation by taking the modulo at each step
def cunr(n, mod=10**6):
    return reduce(lambda a, b: a * b % mod, range(2 * n - 5, 1, -2))


def main(file):
    print(cunr(int(open(file).read())))
