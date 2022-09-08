# Overlap Graphs

from .helpers import Parser
from itertools import permutations


def overlap_graph(seqs, n=3):
    for pair in permutations(seqs, 2):
        if pair[0].seq.endswith(pair[1].seq[:n]):
            yield (pair[0].id, pair[1].id)


def main(file):
    fa = Parser(file).fastas()
    for i in overlap_graph(fa):
        print(*i)
