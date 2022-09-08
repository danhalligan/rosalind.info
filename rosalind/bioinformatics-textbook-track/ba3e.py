# Construct the De Bruijn Graph of a String

from .ba3d import dbru


def main(file):
    seqs = open(file).read().splitlines()
    for k, v in dbru(seqs).items():
        print(k + " -> " + ",".join(v))
