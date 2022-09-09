# The Founder Effect and Genetic Drift

from math import log10
from .helpers import Parser
from .wfmd import wf_model


def foun(n, m, a):
    """The Founder Effect and Genetic Drift"""
    for g in range(1, m + 1):
        yield [log10(wf_model(2 * n, i, g)[0]) for i in a]


def main(file):
    """The Founder Effect and Genetic Drift"""
    l1, l2 = Parser(file).lines()
    n, m = [int(x) for x in l1.split()]
    a = [int(x) for x in l2.split()]
    for x in foun(n, m, a):
        print(*[f"{f:.12f}" for f in x])
