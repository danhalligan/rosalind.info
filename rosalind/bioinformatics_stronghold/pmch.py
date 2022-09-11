# Perfect Matchings and RNA Secondary Structures

from math import factorial, prod
from .helpers import Parser


def pmch(seq):
    return prod([factorial(seq.count(x)) for x in "AG"])


def main(file):
    print(pmch(Parser(file).seqs()[0]))
