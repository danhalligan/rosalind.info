# Introduction to Random Strings

from math import log10
from .helpers import Parser


def logp(x, gc):
    """log10 probability of observing x given GC content"""
    return log10({"G": gc, "C": gc, "A": 1 - gc, "T": 1 - gc}[x] / 2)


def prob(seq, arr):
    """Introduction to Random Strings"""
    return [sum([logp(x, gc) for x in seq]) for gc in arr]


def main(file):
    seq, arr = Parser(file).lines()
    arr = [float(x) for x in arr.split()]
    res = [f"{x:.3f}" for x in prob(seq, arr)]
    print(*res)
