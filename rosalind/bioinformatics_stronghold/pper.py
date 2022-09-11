# Partial Permutations

from functools import reduce
from .helpers import Parser


def main(file):
    n, k = Parser(file).ints()
    print(reduce(lambda p, i: (p * i) % 1_000_000, range(n, n - k, -1)))
