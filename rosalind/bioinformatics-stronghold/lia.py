# Independent Alleles

from .helpers import Parser
from math import comb


def dbinom(x, size, prob):
    """Binomial density"""
    return comb(size, x) * prob**x * (1.0 - prob) ** (size - x)


def pbinom(q, size, prob):
    """Binomial distribution function"""
    return sum([dbinom(x, size, prob) for x in range(0, q + 1)])


def lia(k, n):
    return 1 - pbinom(n - 1, 2**k, 0.25)


def main(file):
    print(round(lia(*Parser(file).ints()), 3))
