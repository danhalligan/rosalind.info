from .ba5d import topological_order
from collections import defaultdict


def parse_graph(graph):
    g = defaultdict(list)
    for edge in graph:
        x, nodes = edge.split(" -> ")
        for y in nodes.split(","):
            g[x] += [{"n": y, "w": 0}]
    return g


def main(file):
    graph = open(file).read().splitlines()
    print(*sorted(topological_order(parse_graph(graph))), sep=", ")
