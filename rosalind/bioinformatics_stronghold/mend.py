# Inferring Genotype from a Pedigree

from .nwck import parse_newick
import numpy as np


# For two vectors of parent genotype probabilities, compute child vector
def child_gt(a, b):
    cross = np.array(
        [
            [[1, 0, 0], [0.5, 0.5, 0], [0, 1, 0]],
            [[0.5, 0.5, 0], [0.25, 0.5, 0.25], [0, 0.5, 0.5]],
            [[0, 1, 0], [0, 0.5, 0.5], [0, 0, 1]],
        ]
    )
    return np.dot(a, np.dot(b, cross))


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
