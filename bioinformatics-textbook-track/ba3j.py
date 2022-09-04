# Reconstruct a String from its Paired Composition

from collections import defaultdict
from .ba3b import genome_path
from .ba3g import eulerian_path


def dbru_paired(pairs):
    g = defaultdict(list)
    for x in pairs:
        p = tuple([x[0][:-1], x[1][:-1]])
        s = tuple([x[0][1:], x[1][1:]])
        g[p].append(s)
    return g


def string_from_paired_composition(pairs, k, d):
    path = eulerian_path(dbru_paired(pairs))
    a = genome_path([x[0] for x in path])
    b = genome_path([x[1] for x in path])
    return a + b[-(k + d) :]


def main(file):
    ints, *pairs = open(file).read().splitlines()
    k, d = map(int, ints.split())
    pairs = [x.split("|") for x in pairs]
    print(string_from_paired_composition(pairs, k, d))
