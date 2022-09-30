# Inferring Genotype from a Pedigree

from .nwck import parse_newick
import numpy as np
from itertools import product


def cross(a, b):
    x = {
        (0, 0): [4, 0, 0],
        (0, 1): [2, 2, 0],
        (1, 0): [2, 2, 0],
        (1, 1): [1, 2, 1],
        (0, 2): [0, 4, 0],
        (2, 0): [0, 4, 0],
        (1, 2): [0, 2, 2],
        (2, 1): [0, 2, 2],
        (2, 2): [0, 0, 4],
    }
    return np.array(x[a, b]) / 4


# For two vectors of parent genotype probabilities, compute child vector
def child_gt(a, b):
    return sum(cross(i, j) * a[i] * b[j] for i, j in product([0, 1, 2], repeat=2))


def main(file):
    tree = parse_newick(open(file).read().rstrip())

    # Initial probability vectors based on known genotypes
    gts = {
        "AA": np.array([1, 0, 0]),
        "Aa": np.array([0, 1, 0]),
        "aa": np.array([0, 0, 1]),
    }

    # Compute vector for nodes not in `gts` dict till we've done all
    q = list(tree.keys())
    while q:
        for n in q:
            if all(c["n"] in gts for c in tree[n]):
                q.remove(n)
                gts[n] = child_gt(*[gts[c["n"]] for c in tree[n]])

    print(*[round(x, 3) for x in gts["0"]])
