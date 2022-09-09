# k-Mer Composition

from .helpers import Parser
from itertools import product


def kmer_perm(k):
    set = ["A", "C", "G", "T"]
    perm = list(product(set, repeat=k))
    return ["".join(x) for x in perm]


def main(file):
    """k-Mer Composition"""
    seq = Parser(file).fastas()[0].seq

    # initialise hash with all possible 4-mers permutations
    d = {k: 0 for k in kmer_perm(4)}

    # Run through 4-mer slices of sequence and increment dictionary keys
    for i in range(len(seq) - 3):
        d[seq[i : (i + 4)]] += 1

    print(*d.values())
