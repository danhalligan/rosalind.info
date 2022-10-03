# General Sink

from .helpers import parse_graphs
from .bfs import bfs


def gs(graph):
    for node in graph:
        dist = bfs(graph, node)
        valid = [x >= 0 for x in dist]
        if all(valid):
            return node
    return -1


def main(file):
    graphs = parse_graphs(open(file), directed=True)
    print(*[gs(g) for g in graphs])
