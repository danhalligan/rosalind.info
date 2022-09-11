# Maximum Matchings and RNA Secondary Structures

from math import factorial
from .helpers import Parser


def nPr(n, k):
    """Returns the number of k-permutations of n."""
    return factorial(n) // factorial(n - k)


def mmch(seq):
    """Maximum Matchings and RNA Secondary Structures"""
    au = [seq.count(x) for x in "AU"]
    gc = [seq.count(x) for x in "GC"]
    return nPr(max(au), min(au)) * nPr(max(gc), min(gc))


def main(file):
    print(mmch(Parser(file).seqs()[0]))
