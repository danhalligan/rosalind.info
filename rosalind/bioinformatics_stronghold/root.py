# Counting Rooted Binary Trees
from functools import reduce


# The answer here is just (2n -3)!!
# We can help the calculation by taking the modulo at each step
def root(n, mod=10**6):
    return reduce(lambda a, b: a * b % mod, range(2 * n - 3, 1, -2))


def main(file):
    print(root(int(open(file).read())))
