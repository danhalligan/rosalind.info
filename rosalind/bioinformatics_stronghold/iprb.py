# Mendel's First Law

from .helpers import Parser
from math import comb


def iprb(k, m, n):
    tot = comb(k + m + n, 2)
    poss = comb(k, 2) + k * m + k * n + m * n / 2 + comb(m, 2) * 3 / 4
    return poss / tot


def main(file):
    print(round(iprb(*Parser(file).ints()), 5))
