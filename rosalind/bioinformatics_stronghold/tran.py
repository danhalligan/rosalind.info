# Transitions and Transversions

from .helpers import Parser


def ts(x, y):
    """Transition?"""
    ts = {"A": "G", "C": "T", "G": "A", "T": "C"}
    return x == ts[y]


def tran(seqs):
    """Transitions and Transversions"""
    mm, tr = 0, 0
    for x, y in zip(seqs[0], seqs[1]):
        if x != y:
            mm += 1
            tr += int(ts(x, y))

    return tr / (mm - tr)


def main(file):
    print(round(tran(Parser(file).seqs()), 11))
