# Independent Segregation of Chromosomes

from math import log10
from .helpers import Parser, pbinom


def indc(n):
    """Independent Segregation of Chromosomes"""
    res = [pbinom(2 * n - x, 2 * n, 0.5) for x in range(1, 2 * n + 1)]
    return [round(log10(x), 3) for x in res]


def main(file):
    print(*indc(*Parser(file).ints()))
