# Construct the Overlap Graph of a Collection of k-mers

from itertools import permutations


def overlap_graph(dna):
    """
    Naive algorithm that finds overlaps by considering all pairs
    """
    for pair in permutations(dna, 2):
        if pair[0].endswith(pair[1][:-1]):
            yield (pair[0], pair[1])


def main(file):
    dna = open(file).read().splitlines()
    for a, b in overlap_graph(dna):
        print(a + " -> " + b)
