# Generate All Maximal Non-Branching Paths in a Graph

from .ba3f import parse_edges
from .ba3k import maximal_nonbranching_paths


def main(file):
    edges = open(file).read().splitlines()
    for path in sorted(maximal_nonbranching_paths(parse_edges(edges))):
        print(" -> ".join(path))
