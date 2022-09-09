# Catalan Numbers and RNA Secondary Structures

from functools import cache
from .helpers import Parser


def valid_pair(x, y):
    pair = {"A": "U", "U": "A", "C": "G", "G": "C"}
    return x == pair[y]


@cache
def cat(seq, mod=10**6):
    """Calculate total number of noncrossing perfect matchings"""
    if len(seq) in range(1):
        return 1
    else:
        return (
            sum(
                cat(seq[1:m]) * cat(seq[m + 1 :])
                for m in range(1, len(seq), 2)
                if valid_pair(seq[0], seq[m])
            )
            % mod
        )


def main(file):
    print(cat(*Parser(file).seqs()))
