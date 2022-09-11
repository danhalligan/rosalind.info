# Implement GreedyMotifSearch with Pseudocounts

from .ba2d import greedy_motif_search


def main(file):
    ints, *dna = open(file).read().splitlines()
    k, t = map(int, ints.split())
    print(*greedy_motif_search(dna, k, pc=1), sep="\n")
