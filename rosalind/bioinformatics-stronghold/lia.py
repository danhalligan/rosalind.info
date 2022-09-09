# Independent Alleles

from .helpers import Parser, pbinom


def lia(k, n):
    return 1 - pbinom(n - 1, 2**k, 0.25)


def main(file):
    print(round(lia(*Parser(file).ints()), 3))
