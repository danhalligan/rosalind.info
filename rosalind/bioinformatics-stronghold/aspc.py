# Introduction to Alternative Splicing

from math import comb
from .helpers import Parser


def main(file):
    n, k = Parser(file).ints()
    print(sum([comb(n, x) for x in range(k, n + 1)]) % 1000000)
