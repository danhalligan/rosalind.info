# Phylogeny Comparison with Split Distance

from .ctbl import ctbl
from .nwck import parse_newick


def main(file):
    _, tree1, tree2 = open(file).read().splitlines()
    t1 = set(ctbl(parse_newick(tree1)))
    t2 = set(ctbl(parse_newick(tree2)))
    print(len(t1 ^ t2))
