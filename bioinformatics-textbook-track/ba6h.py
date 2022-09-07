# Implement ColoredEdges

import re
from .ba6f import chromosome2cycle


def ints(x):
    return list(map(int, x.split()))


def colored_edges(P):
    g = {}
    for chromosome in P:
        nodes = chromosome2cycle(chromosome)
        for j in range(len(chromosome)):
            i1 = 2 * j + 1
            i2 = (2 * j + 2) % len(nodes)
            g[nodes[i1]] = nodes[i2]
    return g


def main(file):
    s = open(file).read().rstrip()
    P = [ints(x) for x in re.findall(r"\((.+?)\)", s)]
    edges = [(k, v) for k, v in colored_edges(P).items()]
    print(*edges, sep=", ")
