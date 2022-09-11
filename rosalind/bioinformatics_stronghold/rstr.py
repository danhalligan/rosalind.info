# Matching Random Motifs

from math import exp
from .helpers import Parser


def main(file):
    l1, seq = Parser(file).lines()
    n, x = map(float, l1.split(" "))
    gc = sum([seq.count(x) for x in "GC"])
    lam = ((1 - x) / 2) ** (len(seq) - gc) * (x / 2) ** gc * n
    print(round(1 - exp(-lam), 3))
