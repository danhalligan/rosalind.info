# Construct the De Bruijn Graph of a String

from collections import defaultdict
from .ba1a import substrings


def dbru(seqs):
    g = defaultdict(list)
    for x in seqs:
        g[x[:-1]].append(x[1:])
    return g


def main(file):
    k, seq = open(file).read().splitlines()
    for k, v in dbru(substrings(seq, int(k))).items():
        print(k + " -> " + ",".join(v))
