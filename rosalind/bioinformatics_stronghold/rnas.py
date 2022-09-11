# Wobble Bonding and RNA Secondary Structures

from functools import cache
from .helpers import Parser


def wobble_pair(x, y):
    pair = {"A": "U", "U": "AG", "C": "G", "G": "CU"}
    return x in pair[y]


@cache
def rnas(seq):
    """Wobble Bonding and RNA Secondary Structures"""
    if len(seq) in range(1):
        return 1
    else:
        return rnas(seq[1:]) + sum(
            rnas(seq[1:m]) * rnas(seq[m + 1 :])
            for m in range(4, len(seq))
            if wobble_pair(seq[0], seq[m])
        )


def main(file):
    print(rnas(Parser(file).line()))
