# Counting Disease Carriers

from math import sqrt
from .helpers import Parser


def afrq(a):
    """Counting Disease Carriers"""
    return [2 * sqrt(x) - x for x in a]


def main(file):
    b = afrq(Parser(file).floats())
    print(*[round(x, 3) for x in b])
